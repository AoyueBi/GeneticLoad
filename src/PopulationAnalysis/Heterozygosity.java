package PopulationAnalysis;


import AoUtils.AoFile;
import AoUtils.SplitScript;
import utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

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
        List<String> groupl = new ArrayList<String>(hm.keySet());
        Collections.sort(groupl);
        //***************************** step two  ****************************//

        //程序运行时，输入输出路径设置
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/003_out";
        String taxaDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/001_treeValidatedGroup_byPloidy";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/001";

        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/002_script_SNPbased/sh_heter20200123.sh";

        //java -jar 033_getSNPHeterbySite.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.subset.vcf /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.cultivar.heter.txt /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/Cultivar.txt &;
        //chr1A_vmap2.1.vcf

        try {

            BufferedWriter bw = IOUtils.getTextWriter(scriptS);
            for (int i = 0; i < groupl.size(); i++) {
                String group = groupl.get(i);
                if(group.equals("ExclusionHexaploid") || group.equals("ExclusionTetraploid"))continue;
                if(group.equals("Hexaploid")){
                    String taxafileS = new File(taxaDirS,group+".txt").getAbsolutePath();
                    String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                        String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
                        String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
                        System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS );
                        bw.write("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS);
                        bw.newLine();
                    }
                }
                if(group.equals("Tetraploid")){
                    String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
                    String taxafileS = new File(taxaDirS,group+".txt").getAbsolutePath();
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                        String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
                        String logfileS = new File(logDirS,"log_chr" + chrArr[j] +  "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
                        System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS );
                        bw.write("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS);
                        bw.newLine();
                    }


                }
                if(group.equals("Ae.tauschii")){
                    String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                    String taxafileS = new File(taxaDirS,group+".txt").getAbsolutePath();
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                        String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
                        String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
                        System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS );
                        bw.write("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS);
                        bw.newLine();
                    }
                }
            }

            bw.flush();
            bw.close();

            new SplitScript().splitScript(scriptS,"heter",14,3);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
}


