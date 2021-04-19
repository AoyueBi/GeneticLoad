/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import Annotation.AnnotationCrossover;
import AoUtils.AoFile;
import AoUtils.AoMath;
import AoUtils.CalVCF;
import AoUtils.CountSites;
import GermplasmInfo.GermplasmInfo;
import GermplasmInfo.TaxaDB;
import Plot.AoMap;
import Plot.PCA;
import Plot.Tree;
import PopulationAnalysis.*;
import WheatGeneticLoad.*;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.StatUtils;
import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaByte;
import pgl.infra.utils.PStringUtils;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class GeneticLoadEntrance {

    public GeneticLoadEntrance() {
//        this.firstProcess();
//        this.plot();
//        this.infoDB(); // 种质信息库
//        this.DBdeleterious();

        //ternary plot analysis
//        this.geneExpression();
        this.rebuildVCF();


    }

    public void rebuildVCF(){
//        new RebuildVCF();
//        new ScriptHapscanner2();
//        new AoHeterozygosity();
//        new  FilterVCF2();
//        new GermplasmInfo();
//        new VariantsSum();
//        new SIFT();
//        new DeleteriousCountbyIndi(); //根据数据库进行个体Load计算
//        new TaxaDB();  //taxa 类
//        new AoIntrogression();

//        new ScriptHapscanner2(); // 进行indel的hapscanner
//        new  FilterVCF2();

//        new AoWheatTriads();
//        new XPCLR();
//        new DeleteriousXPCLR2();

//        new Fst();
//        new Pi();
//        new TajimaD();
//        new AoP();
//        new FdVSdel();
//        new AoMap();
        new DBgene();
//        new RefBiasEvaluation();
//        new VMap2Cal();

    }




    public void geneExpression(){
//        new GeneExpressionbywheat();

    }



    public void DBdeleterious(){
//        new DeleteriousBiologyAoyue();
//        new EstSFS();
//        new DeleteriousCountbyPop();
//        new HomoeologGenesAnalysis();
    }

    public void infoDB(){
        new GermplasmInfo();
    }


    public void plot() {
        //new CalVCF()
        //new Circos();
//        new Tree();
//        new PCA();

    }

    public void firstProcess() {
        //new MapMake();
        //new Wheat120cleandataProcessor();  //Jiao
        //new Wheat120bamProcessor(); //Jiao
//        new WheatBamDatabase();

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
//        new ScriptHapscanner2();
//        new VariantsSum();
//        new PopGenParaWheat();
        //new SplitScript();
//        new Script();
//        new FilterVCF();
//        new Bin();
//        new AoMath();
//        new TreePreparation();
//        new CalVCF();
//        new BadMutations();
//        new AoHeterozygosity();
//        new Fst();
//        new Pi();
//        new TajimaD();
//        new AnnotationCrossover();
//        new XPCLR();
//        AoMath.topK();
//        new GOanalysis();
//        new DeleteriousXPCLR();
//        new GermplasmInfo();
//        new DBgene();
//        new GeneExpressionbywheat();
    }


    public static void main(String[] args) {
        //ChrConvertionRule c=new ChrConvertionRule(Paths.get("/Users/Aoyue/Documents/Data/wheat/chrConvertionRule.txt"));
//        ChrConvertionRule c = new ChrConvertionRule(Paths.get("/data4/home/aoyue/vmap2/analysis/000_taxaList/chrConvertionRule.txt"));
//        VCF.mergeVCFtoLineage(args[0], args[1], c);
        System.out.println("Aoyue Repository --- Here is the entrance of GeneticLoad!\n");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");

        new GeneticLoadEntrance();
//        new GeneticLoadEntrance(args);
//        new CountSites().filterSNPtoBi_parallel(args[0], args[1]);
//        new FilterVCF2().filter_singleThread(args[0], args[1]);
//        CountSites.mergeVCFbysubgenome(args[0], args[1]);
//        CountSites.mergeVCFtoABsubgenome(args[0], args[1]);
//        CountSites.mergeVCFtoAB_Dsubgenome(args[0], args[1]);


        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");

    }

    public GeneticLoadEntrance(String[] args){
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
//        new CountSites().filterAlleleMaf(args[0], args[1]);
//        new CountSites().filterSNPtoBi(args[0], args[1]);
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
//        new SIFT().reverseRefAltallelebyExonVCF(args[0], args[1]);

//        String infileS="/Users/Aoyue/Documents/chr002.subgenome.maf0.01.SNP_bi.subset.vcf";
//        String outfileS="/Users/Aoyue/Documents/chr002.subgenome.maf0.01.SNP_bi.cultivar.vcf.txt";
//        String taxaList="/Users/Aoyue/Documents/Cultivar.txt";
//        new CalVCF().getSNPHeter(infileS,outfileS,taxaList);
//        new CalVCF().getSNPHeter(args[0], args[1], args[2]);

//        new CountSites().extractVCF(infileS,outfileS,taxaList);
//        new AoHeterozygosity().mkGenotype(args[0], args[1]);
//        new Bin().calwindowstep_ResidualHeterozygosity(args[0], Integer.parseInt(args[1]),Integer.parseInt(args[2]), args[3]);
//        new Fst().mkFstTable(args[0], args[1]);
//        new CalVCF().reheader(args[0], args[1], args[2]);
//        new CalVCF().extractIDHapPosRefAlt(args[0], args[1]);
//        new CalVCF().extractIDHapPosRefAlt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/002_snp");

//        CalVCF.filterHeterinVCF(args[0], Double.parseDouble(args[1]), args[2]);
//        new XPCLR().calDensity(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]),Integer.parseInt(args[3]), Integer.parseInt(args[4]),args[5]);
//        new XPCLR().getGenotypeXPCLR(args[0], args[1], args[2]);
//        new XPCLR().getGenotypeXPCLR_parallele(args[0], args[1], args[2]);
//        new XPCLR().getGenotypeXPCLR_parallele_tetra(args[0], args[1], args[2]);
//        new XPCLR().getGenotypeXPCLR_parallele_diploid(args[0], args[1], args[2]);


//        CalVCF.extractVCFtable(args[0], args[1], args[2]);
//        CountSites.countSitesinFastCallformat(args[0]);

//        new FilterVCF2().filter2(args[0]);
//        new RebuildVCF().checkErrorFromFastCall(args[0]);
//        new XPCLR().step3_getAlleleCountXPCLR_3(args[0],args[1],args[2],args[3]);

//        CalVCF.filterMAFinVCF(args[0],Double.parseDouble(args[1]),args[2]);

    }

}
