package PopulationAnalysis;


import AoUtils.AoFile;

import java.util.HashMap;

/**
 *
 * @author Aoyue
 */
public class Heterozygosity {
    public Heterozygosity(){
        this.scriptSNPbased();
    }

    public void scriptSNPbased(){

        //***************************** step one : 确定其倍性，根据倍性计算 ****************************//
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        new AoFile().readheader(dbfileS);
        HashMap<String,String> hm = new AoFile().getHashMap(dbfileS,8,10);
        System.out.println(hm.keySet() + "\t");


        //程序运行时，输入输出路径设置
        String infileDirS = "";
        String outfileDirS = "";

        //java -jar 033_getSNPHeterbySite.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.subset.vcf /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.cultivar.heter.txt /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/Cultivar.txt &;
        //chr1A_vmap2.1.vcf

    }
    
}


