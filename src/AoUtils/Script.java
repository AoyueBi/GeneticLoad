/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class Script {

    public Script() {
        //System.out.println("");
//        this.splitBwaScript("/Users/Aoyue/Documents/sh_fillterMiss20191120.sh", "sh_filterMiss", 21, 2);
//        this.universalScript();
        //this.removeBadTaxafromVCF();
//        this.cp();
//        this.bgzip_D();
//        this.bgzip_AB();
//        this.bgzip_ABD();

//        this.script_ABD();
//        this.script_AB();
//        this.script_D();

//        this.script_AB_byRef();
//        this.script_ABD_byRef();
//        this.splitScript();
//
//        this.script_bw();
        this.script_sout();

//        this.window_single();
//

    }

    public void window_single(){

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/011_test_checkBrt1Gene/002_changeGeneticPos/15.txt";
        AoFile.readheader(infileS);
        int chrColumn = 100;
        int posIndex = 2;
        int valueIndex = 2;
        double window = 2000;
        double step = 2000;
        String name = new File(infileS).getName().split(".txt")[0] + "_" + window + "window_" + step + "step.txt.gz";
        String parent = new File(infileS).getParent();
        String outfileS = new File(parent,name).getAbsolutePath();
        new AoWinScan().getwindowDistrbution_general(infileS,chrColumn,posIndex,valueIndex,window,step,outfileS);

    }

    /**
     * 滑窗处理
     */
    public void window_parallel(){

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/001_fst";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/003_fst_windowbyJava";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            AoFile.readheader(infileS);
            int chrColumn = 0;
            int posIndex = 1;
            int valueIndex = 4;
            double window = 2000000;
            double step = 1000000;
            String name = new File(infileS).getName().split(".fst")[0] + "_" + window + "window_" + step + "step.txt.gz";
            String parent = new File(infileS).getParent();
            String outfileS = new File(outfileDirS,name).getAbsolutePath();
            new AoWinScan().getwindowDistrbution_general(infileS,chrColumn,posIndex,valueIndex,window,step,outfileS);
        });
    }


    /**
     * 将标准化的结果进行滑窗处理
     */
    public void window(){

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/Cultivar_VS_Landrace_EU_exonRegion_0.0001_100_500.xpclr.txt.gz";
        AoFile.readheader(infileS);
        int chrColumn = 8;
        int posIndex = 9;
        int valueIndex = 10;
        double window = 100000;
        double step = 50000;
        String name = new File(infileS).getName().split(".txt")[0] + "_" + window + "window_" + step + "step.txt.gz";
        String parent = new File(infileS).getParent();
        String outfileS = new File(parent,name).getAbsolutePath();

        System.out.println("nohup java -jar 051_AoWindowScan.jar " + infileS + " " + chrColumn + " " + posIndex + " " + valueIndex + " " + window + " " + step + " " + outfileS + " &" );
    }


    public void script_sout() {
//        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
//        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

//        int[] darray = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("rm -f chr" + chr + "_vmap2.0.vcf.gz");
//            System.out.println("rm -fr output" + chr);
            System.out.print( "\""  + chrArr[i] + "\"" + " " );
        }

    }

    public void script_bw() {
        String outfileS = "";
//        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
//        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
//        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};

//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};


        try {
            BufferedWriter bw = IOUtils.getTextWriter("");
            for (int i = 0; i < chrArr.length; i++) {
                String chr = chrArr[i];
                bw.write("");
                bw.newLine();
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void splitScript(){
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/A_mkmd5.sh",20,5);
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/AB_mkmd5.sh",20,12);
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/ABD_mkmd5.sh",20,22); // 435/20

        for (int i = 424; i < 434; i++) {
//            System.out.println("mv ABD_0" + i + ".bam " + " ABD_0" + (i+1) + ".bam");
            System.out.println("mv ABD_0" + i + ".bam.bai " + " ABD_0" + (i+1) + ".bam.bai");

        }
    }



//vcftools --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/all/chr.ABsubgenome.maf0.005.bi_subset.vcf.gz --keep /data4/home/aoyue/vmap2/analysis/013_subsetvcf/all/BreadWheat_S424.txt --recode --recode-INFO-all --stdout | bgzip -c -@ 40 > /data4/home/aoyue/vmap2/analysis/013_subsetvcf/all/chrA_Bsubgenome.breadWheat_S424.vcf.gz &


    public void script_ABD_byRef() {
//        for (int i = 1; i < 8; i++) {
//            String[] chr = {i + "A", i + "B",i+"D"};
//            for (int j = 0; j < chr.length; j++) {
//                System.out.println("vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr" + chr[j] + "_vmap2.1.vcf --indv PI205738 --recode --stdout | bgzip -c -@ 4 > chr" + chr[j] + "_vmap2.1_heter_SNPbased_Cultivar.vcf.gz &" );
////                System.out.println("vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr" + chr[j] + "_vmap2.1.vcf --indv PI436228 --recode --stdout | bgzip -c -@ 4 > chr" + chr[j] + "_vmap2.1_heter_SNPbased_Landrace_PI436228.vcf.gz &" );
//
//            }
//        }

        for (int i = 1; i < 8; i++) {
            String[] chr = {i + "A", i + "B"};
            for (int j = 0; j < chr.length; j++) {
                System.out.println("vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr" + chr[j] + "_vmap2.1.vcf --indv PI355466 --recode --stdout | bgzip -c -@ 4 > chr" + chr[j] + "_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466.vcf.gz &" );
//                System.out.println("vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr" + chr[j] + "_vmap2.1.vcf --indv PI466987 --recode --stdout | bgzip -c -@ 4 > chr" + chr[j] + "_vmap2.1_heter_SNPbased_WildEmmer_PI466987.vcf.gz &" );

            }
        }

//        for (int i = 1; i < 8; i++) {
//            String[] chr = {i + "D"};
//            for (int j = 0; j < chr.length; j++) {
//                System.out.println("vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr" + chr[j] + "_vmap2.1.vcf --indv AE1211 --recode --stdout | bgzip -c -@ 4 > chr" + chr[j] + "_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz &" );
//
//            }
//        }

    }

    public void script_AB_byRef() {
//        String[] db = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
//        Arrays.sort(db);
        for (int i = 1; i < 8; i++) {
            String[] chr = {i + "A", i + "B"};
            for (int j = 0; j < chr.length; j++) {
                System.out.println(chr[j] + "vcftools");
            }
        }
    }

    /**
     * 目的：将单线程产生的所有log文本合并成一个文件。
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergelogTxt(String infileDirS, String outfileS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            ///读表头，在第4行
            int header = 3; //需要进行修改，表头所在行的索引；
            int deslinenumber = 4; //需要进行修改， 目标行所在的索引；
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < header; i++) { //
                String temp = br.readLine();
            }
            bw.write(br.readLine()); //第四行是表头
            bw.newLine();

            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = null;
                //int chrint = Integer.parseInt(fs[i].getName().substring(3, 6));
                int cnt = 0;
                String goal = null;
                for (int j = 0; j < 50; j++) {  //每个log文件有7行，我们只要第5行的数据
                    temp = br.readLine();
                    if (j == deslinenumber) { //目标行的索引
                        StringBuilder sb = new StringBuilder();
                        goal = temp;
                        sb.append(goal);
                        bw.write(sb.toString());
                        bw.newLine();

                    }
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + goal); //print the goal lines
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * find -name "*.vcf" | cut -c3- ; 本地获取运行脚本，输出在netbeans的output界面。
     */
    public void universalScript() {
        try {
            String infileS = "/Users/Aoyue/Documents/a.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
//                String chr = temp.substring(3, 6);
//                System.out.println("mv " + temp + " chr" + chr + "_vmap2_subset0.001.vcf.gz");
//                System.out.println("bgzip -c -@ 4 " + temp + " > " + "/data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/001_byPloid/007_refHexaploid_bgzip/" + temp + ".gz &");
//                System.out.println("tabix -p vcf " + temp + ".gz &");
                System.out.println("cp -f /data2/aoyue/vmap2_AB_project/003_mergefile/AB/" + temp + ".rmdup.bam /mnt/usb/wheatVMapII_ABgenome_merge_rmdupBam_001_1-31/ && cp -f /data2/aoyue/vmap2_AB_project/003_mergefile/AB/" + temp + ".rmdup.bam.bai /mnt/usb/wheatVMapII_ABgenome_merge_rmdupBam_001_1-31/");

            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void script_AB() {
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            if (Arrays.binarySearch(db, i) < 0) { //是属于AB的
//                System.out.println("java -jar 025_cntSitesinMergedVCFtoPop.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/010_miss0.2byPop/chr" + chr + "_vmap2.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/EmmerWheat_S187.txt > log_025/emmer/log_cntSitesinMergedVCFtoPop_chr" + chr + "_tetraploid.txt 2>&1 &");
//                System.out.println("java -jar 028_extractVCF.jar  /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/001_byPloid/tetraploid/chr" + chr + "_vmap2.1_tetraploid.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/EmmerWheat_S187.txt > log_028/log_extractVCF_chr" + chr + "_tetraploid20191107.txt 2>&1");
//                System.out.println("java -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/002_bySubspecies/003_WildEmmer/chr" + chr + "_vmap2.1_WildEmmer.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/002_bySubspecies/tetraploid/Wild_emmer.txt > log_028/001_subspecies/log_extractVCF_chr" + chr + "_WildEmmer20191107.txt 2>&1");
//                System.out.println("java -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/002_bySubspecies/004_DomesticatedEmmer/chr" + chr + "_vmap2.1_DomesticatedEmmer.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/002_bySubspecies/tetraploid/Domesticated_emmer.txt > log_028/001_subspecies/log_extractVCF_chr" + chr + "_DomesticatedEmmer20191107.txt 2>&1");
//                System.out.println("java -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/002_bySubspecies/005_FreeThreshingTetraploid/chr" + chr + "_vmap2.1_FreeThreshingTetraploid.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/002_bySubspecies/tetraploid/Free_threshing_tetraploid.txt > log_028/001_subspecies/log_extractVCF_chr" + chr + "_FreeThreshingTetraploid20191107.txt 2>&1");
//                this.splitBwaScript("/Users/Aoyue/Documents/sh_extractVCFbySubspecies20191108.sh", "extractVCFbySubspecies", 15, 11);
            }
        }
    }

    public void script_D() {
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            if (Arrays.binarySearch(db, i) > -1) { //是属于D的
//                System.out.println("java -jar 025_cntSitesinMergedVCFtoPop.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/010_miss0.2byPop/chr" + chr + "_vmap2.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/Ae.tauschii_S36.txt > log_025/tauschii/log_cntSitesinMergedVCFtoPop_chr" + chr + "_diploid.txt 2>&1 &");
//                System.out.println("mv chr" + chr + ".subgenome.vcf ../004_rawMergedVCF_removeBadTaxa_Dsubgenome_threshold1/");
//                System.out.println("java -jar 028_extractVCF.jar  /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/001_byPloid/diploid/chr" + chr + "_vmap2.1_diploid.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/Ae.tauschii_S36.txt > log_028/log_extractVCF_chr" + chr + "_diploid20191107.txt 2>&1");
//                this.splitBwaScript("/Users/Aoyue/Documents/extractVCF20191107.sh", "extractVCF", 18, 3);
                System.out.println("bgzip -@ 6 chr" + chr + ".vcf &" );
            }
        }
    }

    public void script_ABD(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/010_step1";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/017_toLipengChrPosMajorMinor";
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            String infileS = new File(infileDirS,"chr"+chr+"_vmap2.1_AnnoDB.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr"+chr+"_vmap2.1_AnnoDB_majorminor.txt").getAbsolutePath();
//            System.out.println("zcat "+infileS+"|cut -f1,2,5,6 > " + outfileS + " &");
            //zcat chr001_vmap2.1_AnnoDB.txt.gz |head -n 10|cut -f1,2,5,6

            System.out.println("gzip chr" + chr + "_10X_fastcall_hexaploid.vcf &");
        }


    }

    public void script_ABD_file() {
        try {
            String scriptS = "/Users/Aoyue/Documents/AAAAAA.sh"; //将命令写到脚本中，直接可以执行
            BufferedWriter bw = IOUtils.getTextWriter(scriptS);
            String cmd = "";
            for (int i = 1; i < 43; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
//                cmd = "java -jar 029_mkSNPsummary_step1.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/analysis/015_annoDB/010_step1/chr" + chr + "_vmap2.1_AnnoDB.txt.gz";

//                cmd = "java -jar 030_mkSNPsummary_step2.jar /data4/home/aoyue/vmap2/analysis/015_annoDB/010_step1/chr" + chr + "_vmap2.1_AnnoDB.txt.gz /data4/home/aoyue/vmap2/daxing/ancestralAllele/chr" + chr + ".wheat.ancestralAllele.txt /data4/home/aoyue/vmap2/analysis/015_annoDB/011_step2/chr" + chr + "_vmap2.1_AnnoDB_addDAF.txt.gz";
//                cmd = "chr" + chr + "_exon_vmap2.1_reverseRefAlt.vcf";
                cmd = "chr" + chr + "_exon_vmap2.1.vcf";

                //                String cmd = "sudo gunzip –c /data1/publicData/wheat/genotype/VMap/VMapII/VMap2.1/chr" + chr + "_vmap2.1.vcf.gz > /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII_empty/chr" + chr + "_vmap2.1.vcf";
//java -Xms200g -Xmx500g -jar 025_cntSitesinMergedVCFtoPop.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/010_miss0.2byPop/chr036_vmap2.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt > log_025/log_cntSitesinMergedVCFtoPop_chr036.txt 2>&1 &

                //chr001.lineage.vcf
//            System.out.println("java -jar 025_cntSitesinMergedVCFtoPop.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/010_miss0.2byPop/chr" + chr + "_vmap2.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt > log_025/log_cntSitesinMergedVCFtoPop_chr" + chr + ".txt 2>&1 &");
//            System.out.println("mv chr" + chr + ".lineage.vcf chr" + chr + ".subgenome.vcf");
//                System.out.println("mv chr" + chr + "_miss0.2.bi.vcf chr" + chr + "_vmap2.1.vcf");
                // java -Xms50g -Xmx200g -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr001_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/001_byPloid/hexaploid/chr001_vmap2.1_hexaploid.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt > log_028/log_extractVCF_chr001_hexaploid20191107.txt 2>&1 &
//            System.out.println("java -jar 028_extractVCF.jar  /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/001_byPloid/hexaploid/chr" + chr + "_vmap2.1_hexaploid.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt > log_028/log_extractVCF_chr" + chr + "_hexaploid20191107.txt 2>&1 &");
//            System.out.println("java -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/002_bySubspecies/001_Landrace/chr" + chr + "_vmap2.1_Landrace.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/002_bySubspecies/hexaploid/Landrace.txt > /data4/home/aoyue/vmap2/aaPlantGenetics/log_028/001_subspecies/log_extractVCF_chr" + chr + "_Landrace20191107.txt 2>&1");
//            System.out.println("java -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/002_bySubspecies/002_Cultivar/chr" + chr + "_vmap2.1_Cultivar.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/002_bySubspecies/hexaploid/Cultivar.txt > /data4/home/aoyue/vmap2/aaPlantGenetics/log_028/001_subspecies/log_extractVCF_chr" + chr + "_Cultivar20191107.txt 2>&1");
//vcftools --gzvcf /data4/home/aoyue/vmap2/analysis/002_bivcf/d/chr005.Dgenome.bi.vcf.gz --max-missing 0.1 --recode --recode-INFO-all --stdout | bgzip -c -@ 10 > /data4/home/aoyue/vmap2/analysis/003_filterMiss/d/chr005.Dgenome.filterMiss.vcf.gz &
//System.out.println("vcftools --gzvcf chr" + chr + "_vmap2_subset0.001.vcf.gz --max-missing 0.2 --recode --recode-INFO-all --stdout | bgzip -c -@ 10 > /data4/home/aoyue/vmap2/analysis/020_subsetvcf/001_fromMAF0.01byPop/001/chr" + chr + "_vmap2_subset0.001.vcf.gz &");
                bw.write(cmd);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            this.splitBwaScript(scriptS, "mkAnnotation_step2", 9, 5);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void script_ABorD() {
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            if (Arrays.binarySearch(db, i) < 0) { //说明是不属于D的
                System.out.println("");

            } else { //说明是属于D的
                System.out.println("");

            }
        }
    }

    ///******************************************* 分割线 *************************************  
    ///******************************************* 分割线 *************************************  
    ///******************************************* 分割线 *************************************  
    public void bgzip_D() {
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
            if (index > -1) {
                System.out.println("bgzip -@ 10 chr" + chr + ".Dgenome.filtered0.75.vcf");
            }
        }
    }

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
                //System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
                //System.out.println("/data1/programs/bcftools-1.8/bcftools reheader --samples /data4/home/aoyue/vmap2/analysis/017_removeBadTaxa/005_test_reheaderVCF/changeTaxaName.txt --threads 10 /data4/home/aoyue/vmap2/genotype/mergedVCF/005_maf0.01SNP/chr" + chr + ".lineage.maf0.01.SNP.vcf -o /data4/home/aoyue/vmap2/genotype/mergedVCF/006_reheader/chr" + chr + ".lineage.maf0.01.SNP.vcf");
                //System.out.println("rm -f chr" + chr + ".lineage.maf0.01.SNP.vcf");
                //System.out.println("mv chr" + chr + ".lineage.maf0.01.SNP.vcf ../005_maf0.01SNP/");
                //System.out.println("java -jar filterMafbyPop.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf &");
                //System.out.println("java -jar filterMafbyPopHexaTetra.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf");
//System.out.println("java -jar 008_calDepthSDPvalue_singlethread.jar /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/001_fastcall_Dgenome/rawVCF/chr" + chr + ".vcf /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/001_depthDB/chr" + chr + ".Dgenome.depth.txt.gz > log_008/log_calDepthSDPvalue_chr" + chr + "_20191024.txt &");
//                System.out.println("java -jar mergePosList.jar /data4/home/aoyue/vmap2/analysis/011_filterVCF/abd/003_filteredVCF/chr" + chr + ".ABDgenome.filtered0.75.vcf /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/003_filteredVCF/chr" + chr + ".Dgenome.filtered0.75.vcf.gz /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/003_filterVCF_Dgenome/004_mergePos/posAllele/chr" + chr + "_PosAllele.txt.gz > log_mergePosList/log_mergePosList_chr" + chr + "_20191025.txt & ");
//                System.out.println("bgzip -@ 10 chr" + chr + ".ABgenome.filtered0.75.vcf");
                System.out.println("bgzip -@ 20 chr" + chr + ".subgenome.vcf");
            }
        }
    }

    /**
     * chr035.ABDgenome.vcf
     */
    public void bgzip_ABD() {
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
//            System.out.println("bgzip -@ 10 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
//            System.out.println("bgzip -c -@ 10 /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr" + chr + "_vmap2.1.vcf > /data4/home/aoyue/vmap2/genotype/mergedVCF/013_bgzip/chr" + chr + "_vmap2.1.vcf.gz && tabix -p vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_bgzip/chr" + chr + "_vmap2.1.vcf.gz");
//            this.splitBwaScript("/Users/Aoyue/Documents/a.txt", "bgzip_vmap2.1_", 4, 11);
//            System.out.println("bgzip -@ 10 chr" + chr + ".ABDgenome.filtered0.75.vcf");
//            System.out.println("bgzip -@ 10 chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf");
            System.out.println("bgzip -@ 10 chr" + chr + "_miss0.2.vcf");

        }
    }

    public void cp() {
        String[] s = {"LLX", "LGD", "HRV-L1", "HUN-L1", "ITA-C1", "ITA-L1", "MEX-L1"};
        for (int i = 0; i < s.length; i++) {
            System.out.println("cp -Rf /mnt/usb/ABD001/" + s[i] + "_1.fq.gz /data2/sharedData/vmap2/fastq/");
            System.out.println("cp -Rf /mnt/usb/ABD001/" + s[i] + "_2.fq.gz /data2/sharedData/vmap2/fastq/");
        }

    }

    /**
     * find | cut -f2 -d"/" bgzip -c -@ 10 chr005.vcf > chr005.vcf.gz &
     * 如果不写路径的话，会直接压缩覆盖原来的文件;如果写路径，则会重新生成一个文件，原来未压缩的文件依旧存在。前提是bgzip 不加 -c参数
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void bgzip_deprecated(String infileDirS, String outfileDirS, String threads) {
        /**
         * ** need to modify ***
         */
        //===========================
        String cmd = null;
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        Arrays.sort(fs);
        for (int i = 0; i < fs.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append("bgzip -c -@ " + threads + " " + new File(infileDirS, fs[i].getName()).getAbsolutePath() + " > " + new File(outfileDirS, fs[i].getName().replaceFirst(".vcf", ".vcf.gz")).getAbsolutePath() + " &");
            cmd = sb.toString();
            System.out.println(cmd);
        }
    }

    /**
     * 1 2 3 4 5 6 后缀分别是Alineage Blineage Dlineage chr024.Dlineage.vcf
     */
    /**
     * @deprecated
     */
    public void bgzip_lineage_deprecated() {
        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }
        for (int i = 0; i < 42; i++) {
            int j = i + 1;
            String chr = PStringUtils.getNDigitNumber(3, j);
            String lineage = hml.get(j);
            if (lineage.equals("D")) {
                System.out.println("vcftools --gzvcf /data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/chr" + chr + "." + lineage + "lineage.vcf.gz --remove /data4/home/aoyue/vmap2/analysis/017_removeBadTaxa/000_badTaxaList/BadTaxa_Hexa_Diploid_S8.txt --recode --recode-INFO-all --stdout > /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf");

            } else {
                System.out.println("vcftools --gzvcf /data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/chr" + chr + "." + lineage + "lineage.vcf.gz --remove /data4/home/aoyue/vmap2/analysis/017_removeBadTaxa/000_badTaxaList/BadTaxa_Hexa_Tetra_S8.txt --recode --recode-INFO-all --stdout > /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf");

            }
        }
    }

    /**
     *
     * @param infileS
     * @param nameprefix, the script name you wanna
     * @param numfile, the file number you wanna split to
     * @param numcmd, the number of CDM in each file
     * eg:"/Users/Aoyue/Documents/sh_md5_WheatVMapII_ABgenome_fixmatePosBam.sh",
     * "md5_WheatVMapII_ABgenome_fixmateBam_", 20, 32
     */
    public void splitBwaScript(String infileS, String nameprefix, int numfile, int numcmd) {
        //String infileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/bwa_20190705needRERUN.sh";
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/splitScript";
        String parentS = new File(infileS).getParent();
        new File(parentS, "splitScript").mkdirs();
        String outfileDirS = new File(parentS, "splitScript").getAbsolutePath();
        String shfileS = new File(parentS, "sh_split.sh").getAbsolutePath();

        try {
            String[] outfileS = new String[numfile];
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[numfile];
            for (int i = 0; i < outfileS.length; i++) {
                String num = PStringUtils.getNDigitNumber(3, i + 1);
                outfileS[i] = new File(outfileDirS, nameprefix + num + ".sh").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outfileS[i]);
                String temp;
                for (int j = 0; j < numcmd; j++) {
                    if ((temp = br.readLine()) != null) {
                        bw[i].write(temp);
                        bw[i].newLine();
                    }
                }
                bw[i].flush();
                bw[i].close();
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            File[] fs = new File(outfileDirS).listFiles();
            fs = IOUtils.listFilesEndsWith(fs, ".sh");
            Arrays.sort(fs);
            BufferedWriter bw = IOUtils.getTextWriter(shfileS);
            for (int i = 0; i < fs.length; i++) {
                bw.write("sh " + fs[i].getName() + " > log_" + fs[i].getName().split(".sh")[0] + ".txt 2>&1 &");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(String[] args) {
        new Script();

    }

}
