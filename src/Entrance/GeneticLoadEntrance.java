/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import AoUtils.Bin;
import AoUtils.CountSites;
import AoUtils.Script;
import AoUtils.SplitScript;
import Plot.Circos;
import PopulationAnalysis.PopGenPara;

import WheatGeneticLoad.SIFT;
import WheatGeneticLoad.Treetest;
import WheatGeneticLoad.VariantsSum;


import utils.IOUtils;

import java.io.BufferedWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * @author Aoyue
 */
public class GeneticLoadEntrance {

    public GeneticLoadEntrance() {
        this.firstProcess();
        //this.secondProcess();

    }

    public void secondProcess() {
        //new CalVCF();
        //new Circos();
        new Treetest();

    }

    public void firstProcess() {
        //new MapMake();
        //new Wheat120cleandataProcessor();  //Jiao
        //new Wheat120bamProcessor(); //Jiao
        //new WheatBamDatabase();

        /**
         * *************************************
         */
        //new Wheat200cleanDataProcessor(); //Lu200ABD
        //new WheatABandDcleandataProcessor(); //Lu106AB_D
        //new ABDvcfProcessor();
        //new ABvcfProcessor();
        //new DvcfProcessor();
        //new SIFT();
        //new CountSites().countSitesinFastCallformat("/data4/home/aoyue/vmap2/genotype/abd/");
        //new CountSites().mergeTxt("/Users/Aoyue/Documents/d_depth_m", "/Users/Aoyue/Documents/002_chr1D-7D.Dgenome.depth.txt.gz");
        //new FilterVCF().mergePosList("/Users/Aoyue/Documents/chr036.ABDgenome.filtered0.75.vcf.gz", "/Users/Aoyue/Documents/chr036.Dgenome.filtered0.75.vcf.gz", "/Users/Aoyue/Documents/chr036.posAllele.txt");
        //new ScriptHapscanner2();
        //new VariantsSum();
        //new PopGenPara();
        //new SplitScript();

    }

    public static void main(String[] args) {
        System.out.println("Aoyue Repository --- Here is the entrance of GeneticLoad!\n");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        new GeneticLoadEntrance();
        //new VariantsSum().classifySNPs(args[0], args[1]);
        //new VariantsSum().classifySNPs("/Users/Aoyue/Documents/test", "/Users/Aoyue/Documents/out");
        //new VariantsSum().addAncAllele_singlethread(args[0], args[1], args[2]);
        //new VariantsSum().addAncAllele_singlethread("/Users/Aoyue/Documents/chr002.lineage.maf0.005.bi.AnnoDB.addSIFT.txt","/Users/Aoyue/Documents/chr002.wheat.ancestralAllele.txt","/Users/Aoyue/Documents/chr002.addAnc.txt");
       //new Script().bgzip_lineage();
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/001_depth", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/chrall.aveDepth.txt");
        //new Bin().mkBintable2(args[0], args[1], args[2]);
        //new CountSites().filterIndelMaf("/Users/Aoyue/Documents/test", "/Users/Aoyue/Documents/out");
       //new CountSites().filterIndelMaf(args[0], args[1]);
        //new CountSites().extractHapPos(args[0], args[1]);
        //new CountSites().mergesubsetVCF(args[0], args[1]);
        //new CountSites().mergesubsetVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/015_subesetVCF/ab", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/015_subesetVCF/all/chr.ABsubgenome.maf0.005.bi_subset.vcf.gz");
        //new CountSites().mergesubsetVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/015_subesetVCF/ad", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/015_subesetVCF/all/chr.ADsubgenome.maf0.005.bi_subset.vcf.gz");

//new CountSites().calIndiHeter("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/015_subesetVCF/all","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/002_mafHeterMiss");
//new CountSites().calVcfAverageDepth(args[0], args[1]);
        //new CountSites().calSNPHetMissMaf(args[0], args[1]);
//new CountSites().countSitesinFastCallformat(args[0]);
        //new CountSites().filterAlleleMaf(args[0], args[1]);
        //new CountSites().mergeTxt("/Users/Aoyue/Documents/005_bin", "/Users/Aoyue/Documents/snpDensity.txt");
        //new CountSites().subsetVCFRandomParallel(args[0], args[1],args[2]);
        //new FilterVCF().statVcfDepth_SD_PValue("/Users/Aoyue/Documents/ab/", "/Users/Aoyue/Documents/ab_depth/");
        //new FilterVCF().statVcfDepth_SD_PValue(args[0], args[1]);
        //new FilterVCF().statVcfDepth_SD_PValue_singlethread(args[0], args[1]);
        //new Bin().mkBintable2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/008_genomeDistribution/001_binBasedPos", "1000000");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/008_genomeDistribution/001_binBasedPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/008_genomeDistribution/002_mergebinBasedPos/chrAsubgenome.delSNP.changeChrPos.1M.binTable.txt");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/008_genomeDistribution/001_binBasedPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/008_genomeDistribution/002_mergebinBasedPos/chrBsubgenome.delSNP.changeChrPos.1M.binTable.txt");
        //new Script().bgzip_lineage();
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");
        

        
        
        /**
         * ******************************* temporary method
         * *********************************************
         */
        /**
         * ******************************* temporary method
         * *********************************************
         */
        /**
         * ******************************* temporary method
         * *********************************************
         */
        //new SplitScript().splitBwaScript("/Users/Aoyue/Documents/sh_md5_WheatVMapII_ABgenome_fixmatePosBam.sh", "md5_WheatVMapII_ABgenome_fixmateBam_", 20, 32);
        //new GeneticLoadEntrance().test();
        //new CountSites().extractHapPos(args[0], args[1]);
        //new Bin().mkBintable("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/005_mergepos", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/005_bin.1M/", 1000000);
        //new CountSites().changechrPos("/Users/Aoyue/Documents/003_mkHapPos/", "/Users/Aoyue/Documents/004_changeHapPos/");
        //new CountSites().changechrPosOnVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d_changChrPos/");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/004_changeHapPos/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/005_merge/chr_merged.Dgenome.txt.gz");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/011_snpDistribution/d/001_bin/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/011_snpDistribution/d/002_mergeBinfile/chrDgenome.SNPdensity.1Mbwind.txt");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/004_changeHapPos/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/005_mergepos/");
        //new CountSites().mergefile1and2("/Users/Aoyue/Documents/input/", "/Users/Aoyue/Documents/output/");
    }

    public void test() {

        String outfileS = "/Users/Aoyue/Documents/userAdd.txt";

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String[] user = {"aoyue", "changbin", "daxing", "feilu", "guest", "jijin", "jingwang", "junxu", "lipeng", "qianqian", "sharedData", "xiaohan", "xuebo", "yaozhou", "zhiliang"};
            for (int i = 0; i < user.length; i++) {
                bw.write("useradd -d /data1/home/" + user[i] + " -m " + user[i] + "\n"
                        + "passwd " + user[i] + "\n"
                        + "chown " + user[i] + " -R /data1/home/" + user[i] + "\n");

            }
            bw.flush();
            bw.close();

            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

}
