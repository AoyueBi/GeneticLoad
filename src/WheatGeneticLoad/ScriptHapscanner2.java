/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

import AoUtils.*;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.graph.r.Histogram;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class ScriptHapscanner2 {

    public ScriptHapscanner2() {
//        this.mkTaxaRefBam();
        //this.mkParameterchr1_42_ABD();
        //this.mkJavaCmdchr1_42_ABD();

//        this.mkParameterchr1_42_AB();
        //this.mkParameterchr1_42_D();
//        this.mkJavaCmdchr1_42_AB();
        //this.mkJavaCmdchr1_42_D();
        //this.ifDone();
        //new Script().script_local("/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/output/VCF", "10");
//        this.bgzip_AB();
        //this.bgzip_ABD();
//        this.bcftools_merge();

        /**
         * new HapScanner  20200514
         */
//        this.getNewTaxaBamMapAABB();
//        this.checkTaxaAgain();
//        this.getNewTaxaBamMapDD();
//        this.getNewTaxaBamMapAABBDD();
//        this.convertTaxaBamMapformat();
//        this.getHapPos();

//        this.mkParameter();
//        this.mkJavaCmd();
//        this.bgzip();
//        this.bcftools_merge2();

        /**
         * QC
         */
//        this.qualityCheck();
//        this.mergeCheckFile();
//        this.getBinTable();
//        this.addSubspecies();


        /**
         * new Hapscanner 20200817 only for Indel
         */

//        this.mkParameter_forIndel();
//        this.mkJavaCmd_Indel();
//        this.bgzip();
//        this.bcftools_merge2(); // 过滤步骤见 FilterVCF2类


    }


    /**
     * 本方法的目的是：建立28条染色体的java 运行脚本。
     */
    public void mkJavaCmd_Indel() {
        //AABB
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/006_CMD/sh_hapScanner_ABgenome_chr1_42_20200817.sh";
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

        //AABBDD
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/006_CMD/sh_hapScanner_ABDgenome_20200817.sh";
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};

        //DD
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/006_CMD/sh_hapScanner_Dgenome_chr1_42_20200515.sh";
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        String ploidy = AoString.getPloidy(chrArr);
//        ploidy = "abd";

        try {
            //nohup java -Xmx100g -jar TIGER_v1.0.1.jar -a HapScanner -p ./parameters_006_d_hapScanner.txt > log_006_d_hapScanner.txt &
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < chrArr.length; i++) {
                String chr = chrArr[i];
                bw.write("java -Xmx100g -jar TIGER_v1.0.1.jar -a HapScanner -p ./parameters_");
                bw.write(chr);
                bw.write("_" + ploidy + "_hapScanner.txt > log_");
                bw.write(chr);
                bw.write("_" + ploidy + "_hapScanner2.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        // run this script
//        System.out.println("nohup sh sh_hapScanner_Dgenome_chr1_42_20200515.sh > log_hapScanner_Dgenome_chr1_42_20200515.txt 2>&1 &");
        // nohup sh sh_hapScanner_ABDgenome_onlyDsub_20200518.sh > log_hapScanner_ABDgenome_onlyDsub_20200518.txt 2>&1 &
    }


    /**
     * 本方法的目的是：建立28条染色体的parameters文件。
     */
    public void mkParameter_forIndel() {
//        // 参数文件本地输出
//        String outfileDirS = "";
//        // HPC Path set
//        String[] chrArr ={};
//        String taxaBamMapFileS = "";
//        String outDirS = "";

//        // 参数文件本地输出 AABB
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/004_para_ab";
        // HPC Path set
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
        String taxaBamMapFileS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/ab/001_taxaRefBam_ABgenome.txt";
        String outDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/ab/out";

        // 参数文件本地输出 AABBDD
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/003_para_abd";
//        // HPC Path set
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String taxaBamMapFileS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/abd/001_taxaRefBam_ABDgenome.txt";
//        String outDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/abd/out";

        // 参数文件本地输出 DD
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/005_para_d";
//        // HPC Path set
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
//        String taxaBamMapFileS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/d/001_taxaRefBam_Dgenome.txt";
//        String outDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/d/out";

        // no change
        String posAlleleFileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/003_variationLibrary_Indel";
        String posFileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/004_hapPos";
        String ploidy = AoString.getPloidy(chrArr);
        try {
            for (int i = 0; i < chrArr.length; i++) {
                String chr = chrArr[i];
                int chrInt = Integer.parseInt(chr);
                String posAlleleFileS = new File(posAlleleFileDirS,"chr" + chr + "_Indel.posAllele.txt.gz").getAbsolutePath();
                String posFileS = new File(posFileDirS,"chr" + chr + "_Indel.pos.txt.gz").getAbsolutePath();
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_" + ploidy + "_hapScanner.txt").getAbsolutePath();

                BufferedWriter bw = AoFile.writeFile(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("@App:\tHapScanner\n" +
                        "@Author:\tFei Lu\n" +
                        "@Email:\tflu@genetics.ac.cn; dr.lufei@gmail.com\n" +
                        "@Homepage:\thttps://plantgeneticslab.weebly.com/\n" +
                        "\n" +
                        "#HapScanner is used to perform genotyping of diploid species from whole genome sequenceing data, based on an existing genetic variation library.\n" +
                        "#To run and pipeline, the machine should have both Java 8 and samtools installed. The lib directory should stay with TIGER.jar in the same folder.\n" +
                        "#Command line example. java -Xmx100g -jar TIGER.jar -a HapScanner -p parameter_hapscanner.txt > log.txt &\n" +
                        "#To specify options, please edit the the parameters below. Also, please keep the order of parameters.\n" +
                        "\n" +
                        "#Parameter 1: The taxaRefBam file containing information of taxon and its corresponding refernece genome and bam files. The bam file should have .bai file in the same folder\n" +
                        taxaBamMapFileS + "\n" +
                        "\n" +
                        "#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from genetic variation library.\n" +
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C.\n" +
                        "#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n" +
                        posAlleleFileS + "\n" +
                        "\n" +
                        "#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup.\n" +
                        posFileS + "\n" +
                        "\n" +
                        "#Parameter 4: The chromosome which will be scanned.\n" +
                        chrInt + "\n" +
                        "\n" +
                        "#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping are more likely to be genotyped as homozygote when the combined error rate is high.\n" +
                        "0.05\n" +
                        "\n" +
                        "#Parameter 6: The path of samtools\n" +
                        "/data1/programs/samtools-1.8/samtools\n" +
                        "\n" +
                        "#Parameter 7: Number of threads\n" +
                        "32\n" +
                        "\n" +
                        "#Parameter 8: The directory of output\n" +
                        outDirS + "\n");

                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 先进行36号染色体的测试，并将拿到的结果进行indel 和 deletion 拆分，分别画出MAF分布图，看修改代码后的Indel genotype 是否正常。
     *
     */
    public void getIndelandDeletionFile(){
        String infileS = "";
        String indelOutfileS = "";
        String deletionOutfileS = "";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(indelOutfileS);
            BufferedWriter bw2 = AoFile.writeFile(deletionOutfileS);

            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                temp = temp.substring(0, 40); //肯定够
                l = PStringUtils.fastSplit(temp);
                StringBuilder sb = new StringBuilder();
                sb = new StringBuilder(l.get(2));

                l = PStringUtils.fastSplit(temp);
                cnt++;

            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }



    }



    /**
     * 为部分taxa添加亚群信息
     */
    public void addSubspecies(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/002_merge/001_taxa_QC.txt.gz";
//        AoFile.readheader(infileS);
//        String outfileS = "";
//        HashMap<String,String> hm = AoFile.getHashMapStringKey(infileS,0,11);
//        AoFile.addColumbyString(infileS2,0,hm,"Subspecies");


        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/013_indiHeter/merge/heter_indivi.txt";
        AoFile.readheader(infileS);
        String outfileS = "";
        HashMap<String,String> hm = AoFile.getHashMapStringKey(infileS,0,8);
        AoFile.addColumbyString(infileS2,0,hm,"Subspecies");


    }

    public void getBinTable(){

//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.05;
//        double max =0.5;
//        String outfileS = "";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/002_merge/001_site_QC.txt.gz";
//        AoFile.readheader(infileS);

//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.05;
//        double max =0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/003_binTable/maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);

//        int indexGroup = 0;
//        int indexValue = 1;
//        double window = 0.05;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/003_binTable/site_heter.txt";
//        Bin.frequency2_byGroup(infileS,indexGroup,indexValue,window,window,outfileS);

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/002_merge/001_taxa_QC.txt.gz";
        AoFile.readheader(infileS);

        int indexGroup = 3;
        int indexValue = 1;
        double window = 0.01;
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/003_binTable/taxa_heter.txt";
        Bin.frequency2_byGroup(infileS,indexGroup,indexValue,window,window,outfileS);
    }

    public void getCol(){
        String[] in = {"AB","ABD","D"};
        AoColor.genomeType(in);
    }

    /**
     * 将生成的以site和taxa为单位的质控的结果进行合并，使六倍体，四倍体，二倍体在一个文件中
     */
    public void mergeCheckFile(){

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/001";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/002_merge";

        // change every time
//        String suffix = "_site_QC.txt.gz";
        String suffix = "_taxa_QC.txt.gz";

        String outfileS = new File(outfileDirS,"001" + suffix).getAbsolutePath();
        File[] fs = IOUtils.listRecursiveFiles(new File(infileDirS));
        fs = IOUtils.listFilesEndsWith(fs,suffix);
        AoFile.mergeTxt_byFileArray(fs,outfileS);
    }


    /**
     * 对HapScanner 的结果进行质控
     */
    public void qualityCheck () {
        //missing, maf, heterozygous proportion
        String qcDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/001_QC_test/001_qualityCheck";

        String abInVCFDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/001_QC_test/ab";
        String abdInVCFDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/001_QC_test/abd";
        String dInVCFDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/001_QC_test/d";

        this.checkQuality(abInVCFDirS, new File(qcDirS, "ab").getAbsolutePath(), "AB");
        this.checkQuality(abdInVCFDirS, new File(qcDirS, "abd").getAbsolutePath(), "ABD");
        this.checkQuality(dInVCFDirS, new File(qcDirS, "d").getAbsolutePath(), "D");

    }

    private void checkQuality (String vcfDirS, String outDirS, String genomeType) {
        File outDir = new File (outDirS);
        outDir.mkdir();
        int fN = 2;
        int size = 20000; //抽样数量
        List<File> fList = AoFile.getFileListInDir(vcfDirS);
        Collections.sort(fList);

        GenotypeGrid[] gts = new GenotypeGrid[fN];
        int totalSiteCount = 0;
        for (int i = 0; i < gts.length; i++) { //对genotypeTable 进行初始化
            gts[i] = new GenotypeGrid(fList.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
            totalSiteCount+=gts[i].getSiteNumber(); //获取 fN 个文件的总位点数
        }
        GenotypeOperation.mergeGenotypesBySite(gts[0], gts[1]); // 合并两个VCF文件
        GenotypeGrid gt = gts[0];
        TDoubleArrayList missingSite = new TDoubleArrayList(); // calculation 1
        TDoubleArrayList hetSite = new TDoubleArrayList(); // calculation 2
        TDoubleArrayList maf = new TDoubleArrayList(); // calculation 3
        TDoubleArrayList missingTaxon = new TDoubleArrayList(); // calculation 4
        TDoubleArrayList hetTaxon = new TDoubleArrayList(); // calculation 5

        int step = gt.getSiteNumber()/size;
        for (int i = 0; i < gt.getSiteNumber(); i+=step) {
            missingSite.add(((double) gt.getMissingNumberBySite(i)/gt.getTaxaNumber()));
            hetSite.add(gt.getHeterozygousProportionBySite(i));
            maf.add(gt.getMinorAlleleFrequency(i));
        }
        for (int i = 0; i < gt.getTaxaNumber(); i++) {
            missingTaxon.add((double)gt.getMissingNumberByTaxon(i)/gt.getSiteNumber());
            hetTaxon.add(gt.getHeterozygousProportionByTaxon(i));
        }

        String siteQCfileS = new File(outDirS,genomeType + "_site_QC.txt.gz").getAbsolutePath();
        String taxaQCFileS = new File (outDirS, genomeType+"_taxa_QC.txt.gz").getAbsolutePath();

        try {
            BufferedWriter bw = AoFile.writeFile(siteQCfileS);
            bw.write("GenomeType\tHeterozygousProportion\tMissingRate\tMaf");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < missingSite.size(); i++) {
                sb.setLength(0);
                sb.append(genomeType).append("\t").append(hetSite.get(i)).append("\t").append(missingSite.get(i)).append("\t").append(maf.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }


        try {
            BufferedWriter bw = AoFile.writeFile(taxaQCFileS);
            bw.write("Taxa\tHeterozygousProportion\tMissRate\tGenomeType");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < gt.getTaxaNumber(); i++) {
                sb.setLength(0);
                sb.append(gt.getTaxonName(i)).append("\t").append(hetTaxon.get(i)).append("\t").append(missingTaxon.get(i)).append("\t").append(genomeType);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        //        String siteMissingPdf = new File (outDirS, genomeType+"_site_missing.pdf").getAbsolutePath();
//        String siteHetPdf = new File (outDirS, genomeType+"_site_het.pdf").getAbsolutePath();
//        String taxaMissingPdf = new File (outDirS, genomeType+"_taxa_missing.pdf").getAbsolutePath();
//        String taxaHetPdf = new File (outDirS, genomeType+"_taxa_het.pdf").getAbsolutePath();
//        String mafPdf = new File (outDirS, genomeType+"__site_maf.pdf").getAbsolutePath();
//        String taxaHetFileS = new File (outDirS, genomeType+"_site_het.txt").getAbsolutePath();

//        try {
//            BufferedWriter bw = IOUtils.getTextWriter(taxaHetFileS);
//            bw.write("Taxa\tHeterozygousProportion");
//            bw.newLine();
//            StringBuilder sb = new StringBuilder();
//            for (int i = 0; i < gt.getTaxaNumber(); i++) {
//                sb.setLength(0);
//                sb.append(gt.getTaxonName(i)).append("\t").append(hetTaxon.get(i));
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            bw.flush();
//            bw.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }

//        Histogram h = new Histogram(missingSite.toArray());
//        h.setTitle(genomeType+"_missing"+"_bySite");
//        h.setXLab("Missing rate");
//        h.setYLab("Proportion");
//        h.saveGraph(siteMissingPdf);
//        h = new Histogram(hetSite.toArray());
//        h.setTitle(genomeType+"_heterozygous"+"_bySite");
//        h.setXLab("Heterozygous proportion by site");
//        h.setYLab("Proportion");
//        h.saveGraph(siteHetPdf);
//        h = new Histogram(maf.toArray());
//        h.setTitle(genomeType+"_MAF"+"_bySite");
//        h.setXLab("MAF");
//        h.setYLab("Proportion");
//        h.saveGraph(mafPdf);
//        h = new Histogram(missingTaxon.toArray());
//        h.setTitle(genomeType+"_missing"+"_byTaxa");
//        h.setXLab("Missing rate");
//        h.setYLab("Proportion");
//        h.saveGraph(taxaMissingPdf);
//        h = new Histogram(hetTaxon.toArray());
//        h.setTitle(genomeType+"_heterozygous"+"_byTaxa");
//        h.setXLab("Heterozygous proportion by taxa");
//        h.setYLab("Proportion");
//        h.saveGraph(taxaHetPdf);
    }

    /**
     * 本方法的目的是进行ABD AB D VCF文件的合并,写成脚本形式 bcftools merge -m all --force-samples
     * -f PASS,.
     * /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/ff/out_hapscanner/VCF/chr001_2.vcf.gz
     * /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/ff/out_hapscanner/VCF/chr001.vcf.gz
     * -o
     * /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/ff/out_hapscanner/VCF/merge.vcf
     */
    public void bcftools_merge2() {

//        String abdFileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/abd/out/VCF";
//        String abFileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/ab/out/VCF";
//        String dFileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/d/out/VCF";
//        String mergedFileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF";

        String abdFileDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/abd/out/VCF";
        String abFileDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/ab/out/VCF";
        String dFileDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/d/out/VCF";
        String mergedFileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/106_MergedIndel_fromHapScanner";

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < 42; i++) {
            String chr = chrArr[i];
            String ABorD = AoString.ABorD(chr);
            String abdPath = new File(abdFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String path2 = null;
            if (ABorD.equals("ab")){
                path2 = new File(abFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            }else {
                path2 = new File(dFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            }

            String mPath = new File(mergedFileDirS, "chr" + chr + ".vcf").getAbsolutePath();
            //线程设置 10 8
            System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 8 " + abdPath + " " + path2 + " -o " + mPath + " &");
        }
    }


    public void bgzip() {

//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
            System.out.println("bgzip -@ 4 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
        }
    }


    /**
     * 本方法的目的是：建立28条染色体的java 运行脚本。
     */
    public void mkJavaCmd() {
        //AABB
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/006_CMD/sh_hapScanner_ABgenome_chr1_42_20200515.sh";
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

        //AABBDD
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/006_CMD/sh_hapScanner_ABDgenome_chr1_42_20200531.sh";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/006_CMD/sh_hapScanner_ABDgenome_onlyDsub_20200518.sh";

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};

        //DD
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/006_CMD/sh_hapScanner_Dgenome_chr1_42_20200515.sh";
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        String ploidy = AoString.getPloidy(chrArr);
        ploidy = "abd";

        try {
            //nohup java -Xmx100g -jar TIGER_v1.0.1.jar -a HapScanner -p ./parameters_006_d_hapScanner.txt > log_006_d_hapScanner.txt &
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < chrArr.length; i++) {
                String chr = chrArr[i];
                bw.write("java -Xmx100g -jar TIGER_v1.0.1.jar -a HapScanner -p ./parameters_");
                bw.write(chr);
                bw.write("_" + ploidy + "_hapScanner.txt > log_");
                bw.write(chr);
                bw.write("_" + ploidy + "_hapScanner2.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        // run this script
//        System.out.println("nohup sh sh_hapScanner_Dgenome_chr1_42_20200515.sh > log_hapScanner_Dgenome_chr1_42_20200515.txt 2>&1 &");
        // nohup sh sh_hapScanner_ABDgenome_onlyDsub_20200518.sh > log_hapScanner_ABDgenome_onlyDsub_20200518.txt 2>&1 &
    }


    /**
     * 本方法的目的是：建立28条染色体的parameters文件。
     */
    public void mkParameter() {
//        // 参数文件本地输出
//        String outfileDirS = "";
//        // HPC Path set
//        String[] chrArr ={};
//        String taxaBamMapFileS = "";
//        String outDirS = "";

//        // 参数文件本地输出 AABB
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/004_para_ab";
//        // HPC Path set
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String taxaBamMapFileS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/ab/001_taxaRefBam_ABgenome.txt";
//        String outDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/ab/out";

//        // 参数文件本地输出 AABBDD
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/003_para_abd";
//        // HPC Path set
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String taxaBamMapFileS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/abd/001_taxaRefBam_ABDgenome.txt";
//        String outDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/abd/out";

        // 参数文件本地输出 DD
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/005_para_d";
        // HPC Path set
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
        String taxaBamMapFileS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/d/001_taxaRefBam_Dgenome.txt";
        String outDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/d/out";

        // no change
        String posAlleleFileDirS = "/data4/home/aoyue/vmap2/feilu/variationLibrary/";
        String posFileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/pileup_library/hapPos";
        String ploidy = AoString.getPloidy(chrArr);
        try {
            for (int i = 0; i < chrArr.length; i++) {
                String chr = chrArr[i];
                int chrInt = Integer.parseInt(chr);
                String posAlleleFileS = new File(posAlleleFileDirS,"chr" + chr + ".lib.txt.gz").getAbsolutePath();
                String posFileS = new File(posFileDirS,"chr" + chr + ".lib_HapPos.txt.gz").getAbsolutePath();
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_" + ploidy + "_hapScanner.txt").getAbsolutePath();

                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("@App:\tHapScanner\n" +
                        "@Author:\tFei Lu\n" +
                        "@Email:\tflu@genetics.ac.cn; dr.lufei@gmail.com\n" +
                        "@Homepage:\thttps://plantgeneticslab.weebly.com/\n" +
                        "\n" +
                        "#HapScanner is used to perform genotyping of diploid species from whole genome sequenceing data, based on an existing genetic variation library.\n" +
                        "#To run and pipeline, the machine should have both Java 8 and samtools installed. The lib directory should stay with TIGER.jar in the same folder.\n" +
                        "#Command line example. java -Xmx100g -jar TIGER.jar -a HapScanner -p parameter_hapscanner.txt > log.txt &\n" +
                        "#To specify options, please edit the the parameters below. Also, please keep the order of parameters.\n" +
                        "\n" +
                        "#Parameter 1: The taxaRefBam file containing information of taxon and its corresponding refernece genome and bam files. The bam file should have .bai file in the same folder\n" +
                        taxaBamMapFileS + "\n" +
                        "\n" +
                        "#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from genetic variation library.\n" +
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C.\n" +
                        "#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n" +
                        posAlleleFileS + "\n" +
                        "\n" +
                        "#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup.\n" +
                        posFileS + "\n" +
                        "\n" +
                        "#Parameter 4: The chromosome which will be scanned.\n" +
                        chrInt + "\n" +
                        "\n" +
                        "#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping are more likely to be genotyped as homozygote when the combined error rate is high.\n" +
                        "0.05\n" +
                        "\n" +
                        "#Parameter 6: The path of samtools\n" +
                        "/data1/programs/samtools-1.8/samtools\n" +
                        "\n" +
                        "#Parameter 7: Number of threads\n" +
                        "32\n" +
                        "\n" +
                        "#Parameter 8: The directory of output\n" +
                        outDirS + "\n");

                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }




    public void getHapPos(){
        String infileDirS = "/data4/home/aoyue/vmap2/feilu/variationLibrary";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/pileup_library/hapPos";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/pileup_library/d_hapPos";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
        Arrays.sort(chrArr);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String chr = f.getName().substring(3,6);
                int index = Arrays.binarySearch(chrArr,chr);
                if (index > -1){
                    String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_HapPos.txt.gz").getAbsolutePath();
                    BufferedReader br = AoFile.readFile(infileS);
                    BufferedWriter bw = AoFile.writeFile(outfileS);
                    String header = br.readLine();
                    String temp = null;
                    List<String> l = new ArrayList<>();
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp);
                        bw.write(l.get(0) + "\t" + l.get(1));
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
                    br.close();
                    System.out.println(f.getName() + "\tis completed at " + outfileS);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 第一： 转换 taxaBam文件的格式，从原来一行一个bam文件 --> 一行一个taxa 对应多个bam 文件， bam文件之间用 TAB键隔开
     * 第二： 添加 CS文件，至2个
     */
    public void convertTaxaBamMapformat(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/002_taxaRefBam2";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);

        fsList.parallelStream().forEach(f -> {
            try {
                String refFileS = null;
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName()).getAbsolutePath();
                List<String> taxaList = AoFile.getStringListbySet(infileS,0);
                Collections.sort(taxaList);
                List<String>[] bamList = new List[taxaList.size()];
                for (int i = 0; i < bamList.length; i++) {
                    bamList[i] = new ArrayList<>();
                }
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write("Taxa\tReference\tBams(A list of bams of the taxon, seperated by the delimiter of Tab)");
                bw.newLine();
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String taxa = l.get(0);
                    String bam = l.get(2);
                    refFileS = l.get(1);
                    int index = Collections.binarySearch(taxaList,taxa);
                    bamList[index].add(bam);
                }
                br.close();

                for (int i = 0; i < taxaList.size(); i++) {
                    String taxa = taxaList.get(i);
                    if(taxa.equals("CS")){
                        taxa = "CS-2018";
                    }
                    bw.write(taxa + "\t" + refFileS);
                    for (int j = 0; j < bamList[i].size(); j++) {
                        bw.write("\t" + bamList[i].get(j));
                    }
                    bw.newLine();
                }
                String sub = f.getName().split("_")[2].substring(0,3);
                if (sub.equals("ABD")){
                    bw.write("CS-2017\t" + refFileS + "\t" + "/data3/wgs/bam/ABD/CS_sg_2017_60X.bam");
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getNewTaxaBamMapAABBDD(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/abd/002_taxaRefBam.ABDgenome.manual.addNAFU.txt";
        String bam_Old2NewfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/bam_Old2New.map.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam/001_taxaRefBam_ABDgenome.txt";
        HashMap<String, String> hmold2newBam = AoFile.getHashMapStringKey(bam_Old2NewfileS,0,1);
        String[] removeTaxa = {"IG140057","PI583718","PI534284","Beaqle","Caruton"}; //先去除Taxa,再修改的名字
        Arrays.sort(removeTaxa);
        String refFileS = "/data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa.gz";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header);
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                String oldBam = l.get(2);
                String newBam = hmold2newBam.get(oldBam);
                int index = Arrays.binarySearch(removeTaxa,taxa);
                if (index>-1)continue; //去除坏的taxa
//                if(taxa.equals("CS")){
//                    taxa="CS_mp_2018_8X";
//                }
                bw.write(taxa + "\t" + refFileS + "\t" + newBam);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            AoMath.countCaseInGroup(outfileS,0);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void getNewTaxaBamMapDD(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/d/003_taxaRefBam.Dgenome.addNAFU.txt";
        String bam_Old2NewfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/bam_Old2New.map.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam/001_taxaRefBam_Dgenome.txt";
        HashMap<String, String> hmold2newBam = AoFile.getHashMapStringKey(bam_Old2NewfileS,0,1);
        String[] removeTaxa = {"KU-2071","TA2462","AE430"}; //先去除Taxa,再修改的名字
        Arrays.sort(removeTaxa);
        String refFileS = "/data1/publicData/wheat/reference/v1.0/D/d_iwgscV1.fa.gz";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header);
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                String oldBam = l.get(2);
                String newBam = hmold2newBam.get(oldBam);
                int index = Arrays.binarySearch(removeTaxa,taxa);
                if (index>-1)continue; //去除坏的taxa
                bw.write(taxa + "\t" + refFileS + "\t" + newBam);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            AoMath.countCaseInGroup(outfileS,0);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 再次验证taxaBamMap文件中
     */
    public void checkTaxaAgain(){
//        String infileS ="";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam/ab/001_taxaRefBam_ABgenome.txt";
//        String infileS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam/001_taxaRefBam_Dgenome.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam/001_taxaRefBam_ABDgenome.txt";
        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        List<String> taxaList = AoFile.getStringList(taxaListFileS,0);
        RowTable <String> t = new RowTable<>(infileS);
        int cnt=0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            String taxa = t.getCell(i,0);
            int index = Collections.binarySearch(taxaList,taxa);
            System.out.println(cnt++ + "\t: " + index);
            if (index <0){
                System.out.println(taxa + "\tis not in VMap2.1" );
            }
        }
    }

    /**
     * 结合第一次生成的VMap2.1，将确定的187 Taxa提取出来并进行新的 TaxaBamMap的建立，
     * 第一： bam名字较之前有所变化，第二：Taxa在后期有去除
     */
    public void getNewTaxaBamMapAABB(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/ab/004_taxaRefBam.ABgenome.removeBadTaxa.addNAFU.addS1.txt";
        String bam_Old2NewfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/bam_Old2New.map.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/001_taxaRefBam/ab/001_taxaRefBam_ABgenome.txt";
        HashMap<String, String> hmold2newBam = AoFile.getHashMapStringKey(bam_Old2NewfileS,0,1);
        String[] removeTaxa = {"PI466930","PI466959","PI272522"}; //先去除Taxa,再修改的名字
        Arrays.sort(removeTaxa);
        String refFileS = "/data1/publicData/wheat/reference/v1.0/AB/ab_iwgscV1.fa.gz";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header);
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                String oldBam = l.get(2);
                String newBam = hmold2newBam.get(oldBam);
                int index = Arrays.binarySearch(removeTaxa,taxa);
                if (index>-1)continue; //去除坏的taxa
                if(taxa.equals("PI428082")){ //去除后修改taxa的名字
                    taxa="PI428082_1";
                }
                if (taxa.equals("PI466959_2")){
                    taxa="PI466959";
                }
                bw.write(taxa + "\t" + refFileS + "\t" + newBam);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            AoMath.countCaseInGroup(outfileS,0);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 本方法的目的是进行ABD AB D VCF文件的合并,写成脚本形式 bcftools merge -m all --force-samples
     * -f PASS,.
     * /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/ff/out_hapscanner/VCF/chr001_2.vcf.gz
     * /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/ff/out_hapscanner/VCF/chr001.vcf.gz
     * -o
     * /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/ff/out_hapscanner/VCF/merge.vcf
     */
    public void bcftools_merge() {
//        String abdFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/output/VCF/";
//        String abFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/ab/output/VCF/";
//        String dFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/d/output/VCF/";
//        String mergedFileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/";
        
        String abdFileDirS = "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/004_hapscannerABD/output/VCF/";
        String abFileDirS = "";
        String dFileDirS = "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/005_hapscanner/output/VCF/";
        String mergedFileDirS = "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/006_bcftoolsMerge/";
        
        /**
         * pseudo-code: 1.建立3个lineage的list,然后进行循环，判断：在A lineage下合并，依次类推。
         */
        List<Integer> lA = new ArrayList<>();
        List<Integer> lD = new ArrayList<>();
        //先进行D的建立
        int j = 5;
        lD.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            lD.add(j);
        }
        int k = 6;
        lD.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            lD.add(k);
        }
        //再进行A的建立
        int a = 1;
        lA.add(a);
        for (int i = 0; i < 6; i++) {
            a = a + 6;
            lA.add(a);
        }
        int aa = 2;
        lA.add(aa);
        for (int i = 0; i < 6; i++) {
            aa = aa + 6;
            lA.add(aa);
        }

        Collections.sort(lA);
        Collections.sort(lD);
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            String abdPath = new File(abdFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String abPath = new File(abFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String dPath = new File(dFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            int index = Collections.binarySearch(lD, i);
            int index2 = Collections.binarySearch(lA, i);
            if (index < 0) { //说明是属于AB的
//                if (index2 > -1) { //说明是属于A的
//                    String mPath = new File(mergedFileDirS, "chr" + chr + ".Alineage.vcf").getAbsolutePath();
//                    System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + abPath + " -o " + mPath + " &");
//                } else { //说明是属于B的
//                    String mPath = new File(mergedFileDirS, "chr" + chr + ".Blineage.vcf").getAbsolutePath();
//                    System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + abPath + " -o " + mPath + " &");
//                }
            } else if (index > -1) { //说明是属于D的
                String mPath = new File(mergedFileDirS, "chr" + chr + ".subgenome.vcf").getAbsolutePath();
                System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 10 " + abdPath + " " + dPath + " -o " + mPath + " &");
            }
        }
    }

    public void tabix_ABD() {
        for (int i = 1; i < 13; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            System.out.println("tabix -p vcf chr" + chr + ".vcf.gz &");
        }
    }

    public void bgzip_ABD() {
        for (int i = 1; i < 43; i++) {
            int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(db);
            if (Arrays.binarySearch(db, i) < 0) {
                continue;
            }
            String chr = PStringUtils.getNDigitNumber(3, i);
            System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
        }
    }

    /**
     * 压缩文件bgzip并建立索引
     */
    public void bgzip_AB() {
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(l, i);
            if (index < 0) {
                System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
            }
        }
    }

    /**
     * 本方法的目的是，根据log文件的最后一行文字信息，判断每条染色体是否运行完毕
     */
    public void ifDone() {
        //String infileDirS = "";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/log/abd/";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        try {
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp, " ");
                }
                if (l.contains("completed")) {
                    System.out.println(fs[i].getName() + " is done");
                    //System.out.println("The last line is " + lD.get(0));
                } else {
                    System.out.println(fs[i].getName() + " is not finished");

                }
                br.close();
            }

            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 本方法的目的是：建立14条染色体的java 运行脚本。
     */
    public void mkJavaCmdchr1_42_D() {
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/hapScanner_Dgenome_chr1_42.sh";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/005_hapscanner/hapScanner_Dgenome_chr1_42.sh";
        //nohup java -Xms500g -Xmx500g -jar HapScanner2.jar parameters_001_abd_hapScanner2.txt > log_001_abd_hapScanner2.txt & 
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < 42; i++) {
                int index = Collections.binarySearch(l, i + 1);
                if (index < 0) {
                    continue;
                }
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                bw.write("java -Xms200g -Xmx200g -jar HapScanner2.jar parameters_");
                bw.write(chr);
                bw.write("_d_hapScanner2.txt > log_");
                bw.write(chr);
                bw.write("_d_hapScanner2.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立28条染色体的java 运行脚本。
     */
    public void mkJavaCmdchr1_42_AB() {
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/hapScanner_ABgenome_chr1_42.sh";
        //nohup java -Xms500g -Xmx500g -jar HapScanner2.jar parameters_001_abd_hapScanner2.txt > log_001_abd_hapScanner2.txt & 
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < 42; i++) {
                int index = Collections.binarySearch(l, i + 1);
                if (index > -1) {
                    continue;
                }
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                bw.write("java -Xms500g -Xmx500g -jar HapScanner2.jar parameters_");
                bw.write(chr);
                bw.write("_ab_hapScanner2.txt > log_");
                bw.write(chr);
                bw.write("_ab_hapScanner2.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立14条染色体的parameters文件。
     */
    public void mkParameterchr1_42_D() {
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/para_d_addNAFU/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/005_hapscanner/para_d_addNAFU";
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        try {
            for (int i = 0; i < 42; i++) {
                int index = Collections.binarySearch(l, i + 1);
                if (index < 0) {
                    continue;
                }
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_d_hapScanner2.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("HapScanner2\n"
                        + "Author: Aoyue Bi, Xuebo Zhao, Fei Lu\n"
                        + "Email: biaoyue17@genetics.ac.cn; xuebozhao@genetics.ac.cn; flu@genetics.ac.cn\n"
                        + "Homepage: http://plantgeneticslab.weebly.com/\n"
                        + "#This program is used to genotype whole genome sequenced (WGS) individuals by scanning an existing haplotype library.\n"
                        + "#The usage is java -Xms10g -Xmx20g HapScanner2.jar parameters_hapScanner2.txt > log.txt &\n"
                        + "#Please keep the order of following parameters\n"
                        + "#The taxaRefBam file containing information of taxon and its corresponding reference genome and bam files. The bam file should have .bai file in the same folder. For the situation of one taxon with multiple bams, the bams can be listed by row.\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/012_hapscanner/d/003_taxaRefBam.Dgenome.addNAFU.txt\n"//taxaRefBam路径
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/005_hapscanner/taxaRefBam.Dgenome.addNAFU.txt\n"//taxaRefBam路径
                        + "#The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from haplotype library.\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/004_merge/posAllele/chr" + chr + "_PosAllele.txt.gz\n" //posAllele路径
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/004_mergePos/posAllele/chr" + chr + "_PosAllele.txt.gz\n" //posAllele路径

                        + "#The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup to select pileup sites.\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/004_merge/hapPos/chr" + chr + "_HapPos.txt.gz\n"
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/004_mergePos/hapPos/chr" + chr + "_HapPos.txt.gz\n"
                        + "#The chromosome which will be scanned\n"
                        + (i + 1) + "\n" //染色体号
                        + "#The path of samtools\n"
                        + "/data1/programs/samtools-1.8/samtools\n" //samtools 路径
                        + "#Number of threads\n"
                        + "32\n" //线程大小
                        + "#The directory of output\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/012_hapscanner/d/output/");        
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/005_hapscanner/output/");

                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立28条染色体的parameters文件。
     */
    public void mkParameterchr1_42_AB() {
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/para_ab_addNAFU/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/para_ab_addNAFU/addS1";
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        try {
            for (int i = 0; i < 42; i++) {
                int index = Collections.binarySearch(l, i + 1);
                if (index > -1) {
                    continue;
                }
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_ab_hapScanner2.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("HapScanner2\n"
                        + "Author: Aoyue Bi, Xuebo Zhao, Fei Lu\n"
                        + "Email: biaoyue17@genetics.ac.cn; xuebozhao@genetics.ac.cn; flu@genetics.ac.cn\n"
                        + "Homepage: http://plantgeneticslab.weebly.com/\n"
                        + "#This program is used to genotype whole genome sequenced (WGS) individuals by scanning an existing haplotype library.\n"
                        + "#The usage is java -Xms10g -Xmx20g HapScanner2.jar parameters_hapScanner2.txt > log.txt &\n"
                        + "#Please keep the order of following parameters\n"
                        + "#The taxaRefBam file containing information of taxon and its corresponding reference genome and bam files. The bam file should have .bai file in the same folder. For the situation of one taxon with multiple bams, the bams can be listed by row.\n"
                        + "/data4/home/aoyue/vmap2/analysis/012_hapscanner/ab/004_taxaRefBam.ABgenome.removeBadTaxa.addNAFU.addS1.txt\n"//taxaRefBam路径
                        + "#The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from haplotype library.\n"
                        + "/data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/004_merge/posAllele/chr" + chr + "_PosAllele.txt.gz\n" //posAllele路径
                        + "#The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup to select pileup sites.\n"
                        + "/data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/004_merge/hapPos/chr" + chr + "_HapPos.txt.gz\n"
                        + "#The chromosome which will be scanned\n"
                        + (i + 1) + "\n" //染色体号
                        + "#The path of samtools\n"
                        + "/data1/programs/samtools-1.8/samtools\n" //samtools 路径
                        + "#Number of threads\n"
                        + "16\n" //线程大小
                        + "#The directory of output\n"
                        + "/data4/home/aoyue/vmap2/analysis/012_hapscanner/ab/output/");

                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立42条染色体的java 运行脚本。
     */
    public void mkJavaCmdchr1_42_ABD() {
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/hapScanner_ABDgenome_chr1_42.sh";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/004_hapscannerABD/hapScanner_ABDgenome_chr1_42.sh";
        //nohup java -Xms500g -Xmx500g -jar HapScanner2.jar parameters_001_abd_hapScanner2.txt > log_001_abd_hapScanner2.txt & 
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < 42; i++) {
                int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
                Arrays.sort(db);
                if (Arrays.binarySearch(db, i + 1) < 0) {
                    continue;
                }
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                bw.write("java -Xms500g -Xmx500g -jar HapScanner2.jar parameters_");
                bw.write(chr);
                bw.write("_abd_hapScanner2.txt > log_");
                bw.write(chr);
                bw.write("_abd_hapScanner2.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立42条染色体的parameters文件。
     */
    public void mkParameterchr1_42_ABD() {
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/para_abd/addNAFU/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/004_hapscannerABD/para_abd/";
        try {
            for (int i = 0; i < 42; i++) {
                int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
                Arrays.sort(db);
                if (Arrays.binarySearch(db, i + 1) < 0) {
                    continue;
                }
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_abd_hapScanner2.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("HapScanner2\n"
                        + "Author: Aoyue Bi, Xuebo Zhao, Fei Lu\n"
                        + "Email: biaoyue17@genetics.ac.cn; xuebozhao@genetics.ac.cn; flu@genetics.ac.cn\n"
                        + "Homepage: http://plantgeneticslab.weebly.com/\n"
                        + "#This program is used to genotype whole genome sequenced (WGS) individuals by scanning an existing haplotype library.\n"
                        + "#The usage is java -Xms10g -Xmx20g HapScanner2.jar parameters_hapScanner2.txt > log.txt &\n"
                        + "#Please keep the order of following parameters\n"
                        + "#The taxaRefBam file containing information of taxon and its corresponding reference genome and bam files. The bam file should have .bai file in the same folder. For the situation of one taxon with multiple bams, the bams can be listed by row.\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/002_taxaRefBam.ABDgenome.manual.addNAFU.txt\n"//taxaRefBam路径
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/004_hapscannerABD/002_taxaRefBam.ABDgenome.manual.addNAFU.txt\n"//taxaRefBam路径
                        + "#The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from haplotype library.\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/004_merge/posAllele/chr" + chr + "_PosAllele.txt.gz\n" //posAllele路径
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/004_mergePos/posAllele/chr" + chr + "_PosAllele.txt.gz\n" //posAllele路径        
                        + "#The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup to select pileup sites.\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/004_merge/hapPos/chr" + chr + "_HapPos.txt.gz\n"
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/004_mergePos/hapPos/chr" + chr + "_HapPos.txt.gz\n"
                        + "#The chromosome which will be scanned\n"
                        + (i + 1) + "\n" //染色体号
                        + "#The path of samtools\n"
                        + "/data1/programs/samtools-1.8/samtools\n" //samtools 路径
                        + "#Number of threads\n"
                        + "33\n" //线程大小
                        + "#The directory of output\n"
                        //+ "/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/output/");
                        + "/data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/004_hapscannerABD/output/");

                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkTaxaRefBam() {
        //******************************************* 需要手动选择 *********************************************************************//
        //ABD
//        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_ABD_S373_germplasmInfo.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/taxaRefBam.ABDgenome.txt";

        //AB
//        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_AB_S205_germplasmInfo.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/taxaRefBam.ABgenome.txt";
        //D
        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_D_S35_germplasmInfo.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/001_taxaRefBam.Dgenome.txt";

//******************************************* 需要手动选择 *********************************************************************//
        String abdBamDirS = "/data3/wgs/bam/ABD/";
        String abBamDirS = "/data3/wgs/bam/AB/";
        String dBamDirS = "/data3/wgs/bam/D/";

        String refABD = "/data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa.gz";
        String refAB = "/data1/home/aoyue/wheatRef_v1.0/AB/ab_iwgscV1.fa.gz";
        String refD = "/data1/home/aoyue/wheatRef_v1.0/D/d_iwgscV1.fa.gz";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine(); //read header
            bw.write("Taxa\tReference\tBamPath(bams of the same taxon can be listed by row)");
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                String id = l.get(0);
                String taxa = l.get(4);
                StringBuilder sb = new StringBuilder();
                //******************************************* 需要手动选择 *********************************************************************//
                //ABD
//                sb.append(taxa).append("\t").append(refABD).append("\t").append(abdBamDirS).append(id).append(".rmdup.bam");
                //AB
//                sb.append(taxa).append("\t").append(refAB).append("\t").append(abBamDirS).append(id).append(".rmdup.bam");
//                //D
                sb.append(taxa).append("\t").append(refD).append("\t").append(dBamDirS).append(id).append(".rmdup.bam");

                //******************************************* 需要手动选择 *********************************************************************//
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
