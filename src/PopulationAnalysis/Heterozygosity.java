package PopulationAnalysis;


import AoUtils.AoFile;

import java.util.HashMap;
import java.util.Set;

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
        HashMap<String,String> hm = new AoFile().getHashMap(dbfileS,10,8);
        System.out.println(hm.entrySet());
        Set<String> group = hm.keySet();
        //***************************** step two  ****************************//

        //程序运行时，输入输出路径设置
        String infileDirS = "";
        String outfileDirS = "";

        //java -jar 033_getSNPHeterbySite.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.subset.vcf /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.cultivar.heter.txt /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/Cultivar.txt &;
        //chr1A_vmap2.1.vcf

        for (int i = 0; i < group.size(); i++) {
            if(group.equals("ExclusionHexaploid") || group.equals("ExclusionTetraploid"))continue;
            if(group.equals("AABBDD")){
                String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    System.out.println("");

                }


            }
            if(group.equals("AABB")){
                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};


            }
            if(group.equals("DD")){
                String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D", "7A","7D"};


            }


        }




    }
    
}


