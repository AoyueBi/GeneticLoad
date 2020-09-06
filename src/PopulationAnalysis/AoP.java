package PopulationAnalysis;

import AoUtils.AoFile;
import AoUtils.AoWinScan;
import AoUtils.SplitScript;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;

/**
 * @author AoyueBi
 * @data 2020-09-06 17:19
 */
public class AoP {

    public AoP(){
        this.mkPCommand();

    }

    public void window(){

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/007_pi";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/008_pi/";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            AoFile.readheader(infileS);
            int chrColumn = 0;
            int posIndex = 1;
            int valueIndex = 4;
            double window = 2000000;
            double step = 1000000;
            String name = new File(infileS).getName().split(".pi")[0] + "_" + window + "window_" + step + "step.txt";
            String parent = new File(infileS).getParent();
            String outfileS = new File(outfileDirS,name).getAbsolutePath();

            new AoWinScan().getwindowDistrbution_general(infileS,chrColumn,posIndex,valueIndex,window,step,outfileS);
        });

    }


    /**
     * 计算每个位点在哈温平衡下的P值
     */
    public void mkPCommand(){

        //        //local path
//        String grouplocalDirS = "";
//        //HPC path
//        String infileDirS = "";
//        String groupDirS = "";
//        String outfileDirS = "";
//        //para
//        String window = "";
//        String step = "";
        //        String scriptS = "";


        /**
         * VMap2.0 after May 1st
         */
        //local path
        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/005_P/000_group";
        //HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/105_VMap2.1ByRef";
        String groupDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/004_P/000_group";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/004_P/002_P/001";
        //para
        int numcmd = 20;
        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/005_P/001_script/P_basedSNP_20200906.sh";



        System.out.println("mkdir 001_srcipt_basedSNP");
        System.out.println("mkdir 002_P");

        List<File> fs = IOUtils.getVisibleFileListInDir(grouplocalDirS);
        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < fs.size(); i++) {
                String groupname = fs.get(i).getName().split(".txt")[0];
                if (groupname.equals("Cultivar") || groupname.equals("Landrace")){
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_basedSNP" ).getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --hardy ");
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }

                }
                else if (groupname.equals("Domesticated_emmer") || groupname.equals("Free_threshing_tetraploid") || groupname.equals("Wild_emmer")){
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_basedSNP").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --hardy ");
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }

                }
                else if (groupname.equals("Ae.tauschii")){
                    String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_basedSNP").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --hardy ");
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }
                }

            }
            bw.flush();bw.close();
            SplitScript.splitScript4(scriptS,numcmd); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
