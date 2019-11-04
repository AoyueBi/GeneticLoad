/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.CountSites;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class ScriptHapscanner2 {

    public ScriptHapscanner2() {
        //this.mkTaxaRefBam();
        //this.mkParameterchr1_42_ABD();
        //this.mkJavaCmdchr1_42_ABD();

        //this.mkParameterchr1_42_AB();
        //this.mkParameterchr1_42_D();
        //this.mkJavaCmdchr1_42_AB();
        //this.mkJavaCmdchr1_42_D();
        //this.ifDone();
        //new Script().script_local("/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/output/VCF", "10");
        //this.bgzip_AB();
        //this.bgzip_ABD();
        //this.bcftools_merge();
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
