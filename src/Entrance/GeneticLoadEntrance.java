/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import AoUtils.Bin;
import AoUtils.CountSites;
import WheatGeneticLoad.Sift;
import java.text.SimpleDateFormat;
import java.util.Date;
import utils.PArrayUtils;

/**
 *
 * @author Aoyue
 */
public class GeneticLoadEntrance {

    public GeneticLoadEntrance() {
        //this.firstProcess();

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
        new Sift();

    }

    public static void main(String[] args) {
        System.out.println("Aoyue Repository --- Here is the entrance of GeneticLoad!\n");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        //new GeneticLoadEntrance();
        //new CountSites().extractHapPos(args[0], args[1]);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");

        /**
         * ******************************* temporary method *********************************************
         */
        /**
         * ******************************* temporary method *********************************************
         */
        /**
         * ******************************* temporary method *********************************************
         */
        //new SplitScript().splitBwaScript("/Users/Aoyue/Documents/sh_md5_WheatVMapII_ABgenome_fixmatePosBam.sh", "md5_WheatVMapII_ABgenome_fixmateBam_", 20, 32);
        //new GeneticLoadEntrance().test();
        //new CountSites().extractHapPos(args[0], args[1]);
        new Bin().mkBintable("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/004_changeHapPos/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/005_bin.1M/", 1000000);
        //new CountSites().changechrPos("/Users/Aoyue/Documents/003_mkHapPos/", "/Users/Aoyue/Documents/004_changeHapPos/");
        //new CountSites().changechrPosOnVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d_changChrPos/");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/004_changeHapPos/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/010_changeChrPos/005_merge/chr_merged.Dgenome.txt.gz");
    }

    public void test() {
        
        int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(600000, 30000);
        for(int i=0; i<bound.length; i++){
//            bound[i][0] = i*30000;
//            bound[i][1] = i*30000;
            System.out.println(bound[i][0]);
        }
        int a = 6;
        
            
    }

}
