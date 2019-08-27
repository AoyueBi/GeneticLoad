/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import AoUtils.CountSites;
import WheatGeneticLoad.FilterVCF;
import java.io.BufferedWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class GeneticLoadEntrance {
    
    public GeneticLoadEntrance() {
        this.firstProcess();
        
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
        //new Script();
        //new CountSites();
        new FilterVCF();
    }
    
    public static void main(String[] args) {
        System.out.println("Aoyue Repository --- Here is the entrance of GeneticLoad!\n");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        new GeneticLoadEntrance();
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
                bw.write("useradd -d /data1/home/" + user[i] + " -m " + user[i]+ "\n"
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
