require "minitest/autorun"
require "bio/db/vcf"
require "bio"
require "stringio"

module TestBio; end
module TestBio::TestDb; end

class TestBio::TestDb::TestVcf < MiniTest::Unit::TestCase
# from http://www.1000genomes.org/node/101
EXAMPLE_VCF = %Q{##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
} + 
%Q{#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3
}.gsub(/ +/, "\t")

  def test_sanity
    Bio::FlatFile.new(Bio::Db::Vcf, StringIO.new(EXAMPLE_VCF)).each_entry do |e|
    end
  end

  def test_simple
    e = Bio::FlatFile.new(Bio::Db::Vcf, StringIO.new(EXAMPLE_VCF)).next_entry
    h = Bio::Db::Vcf::VcfHeader.new(EXAMPLE_VCF.split("\n")[0..17] * "\n" + "\n")
    e.reheader!(h)

    assert_equal("20", e.chrom)
    assert_equal(14370, e.pos)
    assert_equal(%w{rs6054257}, e.id)
    assert_equal("G", e.ref)
    assert_equal(%w{A}, e.alt)
    assert_equal(29.0, e.qual)
    assert_equal(%w{PASS}, e.filter)
    
    assert_equal(3, e.info["NS"])
    assert_equal(true, e.info["H2"])
    assert_equal(true, e.info["DB"])
    assert_equal(14, e.info["DP"])
    assert_equal([0.5], e.info["AF"])

    assert_equal(%{GT:GQ:DP:HQ}.split(":"), e.format)
    assert_equal(%w{NA00001 NA00002 NA00003}, e.genotypes.keys)
    assert_equal([51, 51], e.genotypes["NA00001"]["HQ"])
    assert_equal(["G", "G"], e.genotypes["NA00001"]["GT"])
  end

  def test_count
    count = 0
    Bio::FlatFile.new(Bio::Db::Vcf, StringIO.new(EXAMPLE_VCF)).each_entry do |e|
      count += 1
    end

    assert_equal 5, count
  end

  def test_header_parser
    parser = Bio::Db::Vcf::VcfHeaderParser.new

    r = parser.meta.parse(EXAMPLE_VCF.split("\n")[0] + "\n")
    assert_equal "fileformat", r[:key]
    assert_equal "VCFv4.0", r[:value]

    s = EXAMPLE_VCF.split("\n")[8] + "\n"
    r = parser.info.parse(s)
    assert_equal "AA", r[:id]
    assert_equal "String", r[:type]
    assert_equal "Ancestral Allele", r[:description][:string]
    assert_equal "1", r[:number][:int]

    s = EXAMPLE_VCF.split("\n")[11] + "\n"
    r = parser.filter.parse(s)
    assert_equal "q10", r[:id]
    assert_equal "Quality below 10", r[:description][:string]

    s = EXAMPLE_VCF.split("\n")[15] + "\n"
    r = parser.format.parse(s)
    assert_equal "DP", r[:id]
    assert_equal "Read Depth", r[:description][:string]

    s = EXAMPLE_VCF.split("\n")[17] + "\n"
    r = parser.header_line.parse(s)
    assert_equal 3, r[:samples].size
    assert_equal "NA00001", r[:samples][0][:sample]
    assert_equal "NA00002", r[:samples][1][:sample]
    assert_equal "NA00003", r[:samples][2][:sample]
  end

  def test_header_transformer
    parser = Bio::Db::Vcf::VcfHeaderParser.new
    r = parser.parse(EXAMPLE_VCF.split("\n")[0..17] * "\n" + "\n")
    r = Bio::Db::Vcf::VcfHeaderTransformer.new.apply(r)
  end

  def test_header_initialize
    h = Bio::Db::Vcf::VcfHeader.new(EXAMPLE_VCF.split("\n")[0..17] * "\n" + "\n")

    assert_equal 3, h.samples.size
    assert_equal 5, h.meta.size
    assert_equal 6, h.info.size
    assert_equal 4, h.format.size
    assert_equal 2, h.filter.size
  end
end
