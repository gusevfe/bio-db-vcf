require 'bio'
require 'logger'
require 'awesome_print'
require 'parslet'
require 'csv'

$log = Logger.new STDERR

module Bio; end
module Bio::Db; end

def split_quoted(str, sep)
    str.chomp.scan(/(?:"(?:\\.|[^"])*"|[^"#{sep}])+/)
end

class Bio::Db::Vcf
  VERSION = "1.0.0"

  attr_reader :header, :line

  class VcfHeader
    Info = Struct.new(:id, :number, :type, :description)
    Filter = Struct.new(:id, :description)
    Format = Struct.new(:id, :number, :type, :description)

    attr_reader :samples, :meta, :info, :filter, :format

    def initialize(s)
      @meta = Hash.new
      @info = Hash.new
      @filter = Hash.new
      @format = Hash.new
      @samples = []
      parser = VcfHeaderParser.new
      r = parser.parse(s)
      VcfHeaderTransformer.new.apply(r).each do |u|
        type, value = u.first
        if type == :samples
          @samples = value unless value == nil
        elsif type == :meta
          @meta[value[:key]] = value[:value]
        elsif type == :info
          @info[value[:id].to_s] = Info.new(value[:id], value[:number], value[:type], value[:description])
        elsif type == :format
          @format[value[:id].to_s] = Format.new(value[:id], value[:number], value[:type], value[:description])
        elsif type == :filter
          @filter[value[:id].to_s] = Filter.new(value[:id], value[:description])
        end
      end
    end

  end

  class VcfHeaderParser < Parslet::Parser
    rule(:string) { 
      str('"') >> ( str('\\') >> any | str('"').absent? >> any ).repeat.as(:string) >> str('"')
    }

    rule(:integer)    { match('[0-9]').repeat(1).as(:int) }

    rule(:id) {
      match['0-9a-zA-Z_\.'].repeat(1)
    }

    rule(:not_newline) {
      match('[^\n]').repeat(1)
    }

    rule(:newline_or_eof) {
      str("\n") | any.absent?
    }

    rule(:meta) {
      str('##') >> id.as(:key) >> str('=') >> not_newline.as(:value) >> newline_or_eof
    }

    rule(:info_type) {
      str("Integer") |
      str("Float") |
      str("Flag") |
      str("Character") |
      str("String")
    }

    rule(:number) {
      integer | match("[.AG]")
    }

    rule(:info) {
      # ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
      str('##INFO=<ID=') >> id.as(:id) >> str(',Number=') >> number.as(:number) >> str(',Type=') >> info_type.as(:type) >> str(',Description=') >> string.as(:description) >> str(">") >> newline_or_eof
    }

    rule(:filter) {
      # ##FILTER=<ID=ID,Description=”description”>
      str('##FILTER=<ID=') >> id.as(:id) >> str(',Description=') >> string.as(:description) >> str(">") >> newline_or_eof
    }

    rule(:format) {
      # ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
      str('##FORMAT=<ID=') >> id.as(:id) >> str(',Number=') >> number.as(:number) >> str(',Type=') >> info_type.as(:type) >> str(',Description=') >> string.as(:description) >> str(">") >> newline_or_eof
    }

    rule(:header_line) {
      str("#") >> str(%w{CHROM POS     ID        REF    ALT     QUAL FILTER INFO} * "\t") >> (str("\tFORMAT") >> (str("\t") >> id.as(:sample)).repeat(1)).maybe.as(:samples) >> newline_or_eof
    }

    rule(:header) {
       (info.as(:info) | filter.as(:filter) | format.as(:format) | meta.as(:meta)).repeat(1) >> header_line
    }

    root(:header)
  end

  class VcfHeaderTransformer < Parslet::Transform
    rule(:int => simple(:x)) { Integer(x.to_s.sub(/^0*([^0])/, '\1')) }
    rule(:string => simple(:x)) { x }
    rule(:sample => simple(:x)) { x.to_s }
  end

  class VcfSplitter < Bio::FlatFile::Splitter::Default
    def initialize(klass, bstream)
      super(klass, bstream)

      @h = ""

      # Read in the header
      while true
        s = bstream.gets
        break unless s
        if s[0] == '#'
          @h.concat s
        else
          bstream.ungets s
          break
        end
      end

      @h = VcfHeader.new(@h)
    end

    def get_parsed_entry
      super
      parsed_entry.reheader!(@h) if parsed_entry
      parsed_entry
    end
  end

  FLATFILE_SPLITTER = VcfSplitter
  DELIMITER = "\n"

  attr_reader :chrom, :pos, :id, :ref, :alt, :qual, :filter, :info, :format, :genotypes

  def initialize(s)
    @h = nil
    @s = s.split("\t")
 end

  def reheader!(h)
    @h = h

    @chrom = @s[0]
    @pos = Integer(@s[1])
    @id = @s[2].split(";")
    @ref = @s[3]
    @alt = @s[4].split(",")
    @qual = (@s[5] == "." ? "." : Float(@s[5]))
    @filter = @s[6].split(";")
    @info = Hash.new
    split_quoted(@s[7], ";").each do |i|
      i =~ /([0-9a-zA-Z_\.]+)=?(.*)$/
      @info[$1] = format_value(@h.info[$1], $2)
    end

    @genotypes = Hash.new
    unless @h.samples.empty?
      @format = split_quoted(@s[8], ":")

      @h.samples.each_with_index do |sample, i|
        @genotypes[sample] = Hash.new
        split_quoted(@s[9 + i], ":").each_with_index do |v, i|
          f = @format[i]
          v = format_value(@h.format[f], v)
          #
          # Special handling
          if f == "GT"
            v = v.chomp.split(/[|\/]/).map { |u| u == "." ? "." : allele_by_id(Integer(u)) }
          end

          @genotypes[sample][f] = v
        end
      end
    end
  end

  def format_value(u, v)
    if u[:number] == 1
      return format_single_value(u[:type], v)
    elsif u[:number] == 0
      return format_single_value(u[:type], nil)
    else
      split_quoted(v, ",").map do |v|
        format_single_value(u[:type], v)
      end
    end
  end

  def allele_by_id(i)
    i == 0 ? @ref : @alt[i - 1]
  end

  def format_single_value(type, v)
    if type == "Integer"
      if v == "."
        return v
      else
        return Integer(v)
      end
    elsif type == "String" or type == "Character"
      return v
    elsif type == "Float"
      return Float(v)
    elsif type == "Flag"
      return true
    end
  end
end
