/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import AoUtils.Bin;
import Plot.Tree;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * @author Aoyue
 */
public class GeneticLoadEntrance {

    public GeneticLoadEntrance() {
        this.firstProcess();
//        this.secondProcess();

    }



    public void secondProcess() {
        //new CalVCF()
        //new Circos();
        new Tree();

    }

    public void firstProcess() {
        //new MapMake();
        //new Wheat120cleandataProcessor();  //Jiao
        //new Wheat120bamProcessor(); //Jiao
        //new WheatBamDatabase();

        /**
         * *************************************
         */
//        new Wheat200cleanDataProcessor(); //Lu200ABD
        //new WheatABandDcleandataProcessor(); //Lu106AB_D
        //new ABDvcfProcessor();
        //new ABvcfProcessor();
        //new DvcfProcessor();
//        new DataStorage();
//        new SIFT();
//        new CountSites();
        //new ScriptHapscanner2();
//        new VariantsSum();
//        new PopGenParaWheat();
        //new SplitScript();
//        new Script();
//        new FilterVCF();
        new Bin();
//        new AoMath();
//        new TreePreparation();
//        new CalVCF();

//new DeleteriousDB();
    }

    public static void main(String[] args) {
        //ChrConvertionRule c=new ChrConvertionRule(Paths.get("/Users/Aoyue/Documents/Data/wheat/chrConvertionRule.txt"));
//        ChrConvertionRule c = new ChrConvertionRule(Paths.get("/data4/home/aoyue/vmap2/analysis/000_taxaList/chrConvertionRule.txt"));
//        VCF.mergeVCFtoLineage(args[0], args[1], c);
        System.out.println("Aoyue Repository --- Here is the entrance of GeneticLoad!\n");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        new GeneticLoadEntrance();

//        VCF vcf=new VCF("/data4/home/aoyue/vmap2/analysis/020_subsetvcf_fromMAF0.01byPop/002_mergedbySub/chr.lineageA.vcf.gz");
//        vcf.addVCF(new VCF("/data4/home/aoyue/vmap2/analysis/020_subsetvcf_fromMAF0.01byPop/002_mergedbySub/chr.lineageB.vcf.gz"));
//        vcf.write("/data4/home/aoyue/vmap2/analysis/020_subsetvcf_fromMAF0.01byPop/003_all/chr.ABsubgenome.vcf.gz");
        //new VariantsSum().classifySNPs(args[0], args[1]);
        //new VariantsSum().addAncAllele_singlethread(args[0], args[1], args[2]);
        //new Bin().mkBintable2(args[0], args[1], args[2]);
        //new CountSites().filterIndelMaf(args[0], args[1]);
        //new CountSites().extractHapPos(args[0], args[1]);
        //new CountSites().mergesubsetVCF(args[0], args[1]);
        //new CountSites().calVcfAverageDepth(args[0], args[1]);
        //new CountSites().calSNPHetMissMaf(args[0], args[1]);
//        new CountSites().countSitesinFastCallformat(args[0]);
        //new CountSites().filterAlleleMaf(args[0], args[1]);
//        new CountSites().filterAllele(args[0], args[1]);
//        new CountSites().subsetVCF(args[0], args[1],args[2]);
        //new FilterVCF().statVcfDepth_SD_PValue(args[0], args[1]);
        //new FilterVCF().statVcfDepth_SD_PValue_singlethread(args[0], args[1]);
        //new FilterVCF().reINFOHexaTetraPloid(args[0], args[1]);
        //new FilterVCF().reINFOHexaDiPloid(args[0], args[1]);
        //分别测试ABD 和 AB 和 D 用
        //new FilterVCF().filterMafbyPopHexaDi(args[0], args[1]);
        //new FilterVCF().filterMafbyPopHexaTetra(args[0], args[1]);
        //new FilterVCF().statVcfDepth_SD_PValue_singlethread(args[0], args[1]);
        //new FilterVCF().mergePosList(args[0], args[1], args[2]);
        //new FilterVCF().extractHapPosfromposAllele(args[0], args[1]);
//        new FilterVCF().filterMafbyPopHexaDi(args[0], args[1]);
//        new FilterVCF().filterMafbyPopHexaTetra(args[0], args[1]);
//        new FilterVCF().filterMissbyPopHexaDi(args[0], args[1]);
//        new FilterVCF().filterMafbyPopHexaTetra(args[0], args[1]);
//        new CountSites().cntSitesinMergedVCFtoPop(args[0], args[1]);
//        new CountSites().extractVCF(args[0], args[1], args[2]);
//        new VariantsSum().mkSNPsummary_step1(args[0], args[1]);
//        new VariantsSum().mkSNPsummary_step2(args[0], args[1], args[2]);
//        new VariantsSum().getCDSannotation(args[0], args[1]);
//new SplitScript().splitBwaScript(args[0], args[1], Integer.parseInt(args[2]),Integer.parseInt(args[3]));
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");
        /**
         * ******************************* temporary method
         * *********************************************
         */
//        new GeneticLoadEntrance().test();
    }

    public void test() {
        try {
            String sample = "a  b c    d";
            String[] arrays = sample.split(" +");
            for (String s : arrays) {
                System.out.println(s);
            }

//           System.out.println(String.valueOf(cnt) + " SNPs output from ");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
