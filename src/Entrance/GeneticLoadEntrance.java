/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import Annotation.AnnotationCrossover;
import AoUtils.*;
import GeneFetchFromNCBI.Step1GenesFromNCBI;
import GermplasmInfo.GermplasmInfo;
import GermplasmInfo.TaxaDB;
import Plot.AoMap;
import Plot.GOanalysis;
import Plot.PCA;
import Plot.Tree;
import PopulationAnalysis.*;
import WheatGeneticLoad.*;
import WheatVMap2S1000.SampleSize2VariantsDiscovery;
import WheatVMap2S1000.VMap2S1000;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.StatUtils;
import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaByte;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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

        //*** ternary plot analysis ***//
//        this.geneExpression(); //analysis on wheat gene expression and deleterious load
//        this.rebuildVCF(); //new analysis on VMap2.0-2020 (主要)
        this.projectVMap2S1000(); //new analysis on VMap2.0-2021-07

    }

    public void projectVMap2S1000(){
        new VMap2S1000(); //进行VCF的fix, QC, and SNP annotation build
//        new Pi();
//        new Fst();
//        AoFile.readheader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/004_thetaW/002_merge001/angsd_subspecies26_geneRegion.txt.gz");
//        new TajimaD();
//        new SampleSize2VariantsDiscovery();

//        String infileS = "/Users/Aoyue/Documents/chr036_vmap2.1.vcf.gz";
//        String outfileS = "/Users/Aoyue/Documents/chr036_DD.vcf";
//        String taxaList = "/Users/Aoyue/Documents/test/taxalist.txt";
//        CalVCF.extractVCFforANGSDinput(infileS,outfileS,taxaList);

//        int chrint = 1;
//        int posint = 1172094;
//        String chr = RefV1Utils.getChromosome(chrint,posint);
//        int pos = RefV1Utils.getPosOnChromosome(chrint,posint);
//        System.out.println("chr "+ chr);
//        System.out.println("pos " + pos);

//        new Step1GenesFromNCBI();

    }


    public void rebuildVCF(){
//        new RebuildVCF();
//        new ScriptHapscanner2();
//        new AoHeterozygosity();
//        new  FilterVCF2();
//        new GermplasmInfo(); // 种质信息库
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
//        new DBgene();
//        new RefBiasEvaluation();
//        new VMap2Cal(); //sample size vs variants discovery

    }




    public void geneExpression(){
        new GeneExpressionbywheat();

    }



    public void DBdeleterious(){
//        new DeleteriousBiologyAoyue();
//        new EstSFS();
//        new DeleteriousCountbyPop();
//        new HomoeologGenesAnalysis();
    }

    public void infoDB(){
//        new GermplasmInfo();
    }


    public void plot() {
        //new Circos();
//        new Tree();
//        new PCA();

    }

    public void firstProcess() {
        new MapMake();
        new Wheat120cleandataProcessor();  //Jiao
        new Wheat120bamProcessor(); //Jiao
        new WheatBamDatabase();
        new Wheat200cleanDataProcessor(); //Lu200ABD
        new WheatABandDcleandataProcessor(); //Lu106AB_D
        new ABDvcfProcessor();
        new ABvcfProcessor();
        new DvcfProcessor();
        new SIFT();
        new CountSites();
        new ScriptHapscanner2();
        new VariantsSum();
        new PopGenParaWheat();
        new FilterVCF();
        new AoMath();
        new CalVCF();
        new BadMutations();
        new AoHeterozygosity();
        new Fst();
        new Pi();
        new TajimaD();
        new AnnotationCrossover();
        new XPCLR();
        new GOanalysis();
        new DeleteriousXPCLR();
        new GermplasmInfo();
        new DBgene();
        new GeneExpressionbywheat();
    }


    public static void main(String[] args) {

        System.out.println("Aoyue Repository --- Here is the entrance of GeneticLoad!\n");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        new GeneticLoadEntrance();
//        new GeneticLoadEntrance(args);

        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");

        //ChrConvertionRule c=new ChrConvertionRule(Paths.get("/Users/Aoyue/Documents/Data/wheat/chrConvertionRule.txt"));
//        ChrConvertionRule c = new ChrConvertionRule(Paths.get("/data4/home/aoyue/vmap2/analysis/000_taxaList/chrConvertionRule.txt"));
//        VCF.mergeVCFtoLineage(args[0], args[1], c);

    }

    public GeneticLoadEntrance(String[] args){

        //        new CountSites().filterSNPtoBi_parallel(args[0], args[1]);
//        new FilterVCF2().filter_singleThread(args[0], args[1]);
//        CountSites.mergeVCFbysubgenome(args[0], args[1]);
//        CountSites.mergeVCFtoABsubgenome(args[0], args[1]);
//        CountSites.mergeVCFtoAB_Dsubgenome(args[0], args[1]);

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
//        CalVCF.calMAFcountfromPop(args[0],args[1],args[2]);
//        CalVCF.filterVCFbyPos(args[0],args[1],args[2]);
//        new VMap2S1000().checkAltCaseCount(args[0],args[1]);
//        AoFile.subsetTxt_parallel(args[0],args[1],args[2]);
//        AoFile.subsetTxt_singleStream(args[0],args[1],args[2]);
//        CalVCF.filterMAFinVCF_parallel(args[0],Double.parseDouble(args[1]),args[2]);
//        CalVCF.calAAFFromPop(args[0],args[1],args[2]);
//        new SampleSize2VariantsDiscovery().sampleSize2SNPdiscovery(Integer.parseInt(args[0]),Integer.parseInt(args[1]),args[2],args[3],args[4],args[5]);
//        new CalVCF().extractVCFforANGSDinput(args[0],args[1],args[2]);

//        CountSites.merge1_42to1A_7DandChangeChrPos_txt(args[0],args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
//        CountSites.mergeVCF1_42to1A_7DandChangeChrPos(args[0],args[1]);
//        CountSites.mergeVCF1_42to1A_7DandChangeChrPos_bySuffix(args[0],args[1],args[2]);
//        CountSites.merge1_42to1A_7DandChangeChrPos_txt2(args[0],args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
//        AoFile.filterTxtbyValue(args[0],Integer.parseInt(args[1]),args[2],Double.parseDouble(args[3]),args[4]);

//        AoFile.ChangeChrPos(args[0],Integer.parseInt(args[1]),Integer.parseInt(args[2]),args[3]);
//        AoWheatTriads.getTriadsModel(args[0],args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]));
//        CalVCF.getGenoTable_includeChrPosRefAlt(args[0],args[1],args[2]);
//        CalVCF.extractVCFtable(args[0],args[1]);

        new VMap2S1000().getProportion(args[0],args[1],args[2],args[3]);

    }

}
