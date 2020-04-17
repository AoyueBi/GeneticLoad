/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PopulationAnalysis;

import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class PopGenParaWheat {

    public PopGenParaWheat() {
//        this.getTaxaSet();
//        this.mkFstCommandbasedSNP();
//        new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/sh", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/fst_basedSNP.sh");
//        new Script().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/fst_basedSNP.sh", "fst_basedSNP", 30, 7);
//        this.fstTable();
//        this.mkFstCommandbasedwinndow();
//        this.mergeSh("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/002_script_based10kbWindow/001_total", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/002_script_based10kbWindow/002_merge_sh/fst_basedWindow20191126.sh");
//        new Script().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/002_script_based10kbWindow/002_merge_sh/fst_basedWindow20191126.sh", "fst_basedWindow_", 30, 6);

//        this.mkPiCMDbasedWindow();
//        this.mergeSh("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/001_script_based100kbWindow/001_total", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/001_script_based100kbWindow/002_merge_sh/pi_based100kbWindow20191126.sh");
//        new Script().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/001_script_based100kbWindow/002_merge_sh/pi_based100kbWindow20191126.sh", "pi_based100kbWindow_", 5, 18);
//   this.mkTajimaDbasedWindow();
//   this.mergeSh("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/001_script_based100kbWindow/001_total", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/001_script_based100kbWindow/02_merge_sh/tajimaD_based100kbWondow20191127.sh");
//   new Script().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/001_script_based100kbWindow/02_merge_sh/tajimaD_based100kbWondow20191127.sh", "tajimad100kbWindow_", 50, 2);
//    this.mergePopfile();
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/003_mergeDsub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/004_mergAll/tajimaD_all_based100kbWindow.txt");

//        this.mkPiCMDbasedWindow();
//        this.mergeandsplitSh("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/003_script_based100kbWindow_byPloidy/001_total", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/003_script_based100kbWindow_byPloidy/002_merge_sh/pi_based100kbWindow20191211.sh", 6, 7);

//        this.mkFstCommandbasedwinndow();
//        this.mergeandsplitSh("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/004_script_based10KbWindow_byPloidy/001_total", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/004_script_based10KbWindow_byPloidy/002_merge_sh/fst_basedWindow20191211.sh", 7, 3);

        /**
         * 按照倍性，将AABBDD AABB DD的文件合并在一起
         */
        this.mergePopfile();

    }
    
    
    
    /**
     * 
     * 将TajimaD等参数的结果合并起来，画在一张图上
     */
    public void mergePopfile() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/001_tajimaD_based100kbWindow_hexa_diploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/003_mergeDsub/Tajima.D_Dsub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/002_tajimaD_based100kbWindow_hexa_tetraploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/003_mergeDsub/Tajima.D_ABsub_based100kbWindow.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/002_tajimaD/003_mergeDsub/Tajima.D_Bsub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/001_pi_based100kbWindow_hexa_diploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/003_merge_sub/Pi_Dsub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/002_pi_based100kbWindow_hexa_tetraploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/003_merge_sub/Pi_Bsub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/001_fst_based100kbWindow_hexa_diploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/003_merge_sub/Fst_Dsub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/002_fst_based100kbWindow_hexa_tetraploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/003_merge_sub/Fst_Asub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/002_fst_based100kbWindow_hexa_tetraploid";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/003_merge_sub/Fst_Bsub_based100kbWindow.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/004_pi_based100kbWindow_byPloidy";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/005_merge_sub_byPloidy/Pi_Asub_based100kbWindow_byPloid.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/005_merge_sub_byPloidy/Pi_Bsub_based100kbWindow_byPloid.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/002_result/005_merge_sub_byPloidy/Pi_Dsub_based100kbWindow_byPloid.txt.gz";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/004_fst_based100kbWindow_byPloidy";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/005_merge_sub_byPloidy/Fst_Asub_based100kbWindow_byPloidy.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/005_merge_sub_byPloidy/Fst_Bsub_based100kbWindow_byPloidy.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/003_result/005_merge_sub_byPloidy/Fst_Dsub_based100kbWindow_byPloidy.txt.gz";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
//        fs = IOUtils.listFilesEndsWith(fs, suffix);
        fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesContains(fs, "A_");
        Arrays.sort(fs);
        try {

            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(infileS);
            ///读表头
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            bw.write(br.readLine() + "\tGroup"); //read header
            bw.newLine();


            int cnttotal = 0;
            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split("_S")[0];
                br = IOUtils.getTextReader(infileS);
                String temp; //read header
                temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp).append("\t").append(group);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS);
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkTajimaDbasedWindow() {
        //vcftools --gzvcf /data2/aoyue/maizeData/hmp321_agp4/hmp321_agpv4_chr1.vcf.gz --keep /data2/aoyue/popGene/002_parameters/000_groups/China_specific.txt --TajimaD 10000 --out /data2/aoyue/popGene/002_parameters/003_tajimaD/000_tajimaDBased10Kwindow/China_specific/China_specific_hmp321_agpv4_chr1 &
//        //脚本文件路径
//        String ScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/001_script_based100kbWindow/hexaTetra";
//        //本地分组文件路径
//        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandTetra";
//        /*vcftools运行文件路径*/
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/";
//        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/001_hexaandTetra";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/003_tajimaD/002_tajimaD_based100kbWindow_hexa_tetraploid";
        
        //脚本文件路径
        String ScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/003_tajimaD/001_script_based100kbWindow/hexaDi";
        //本地分组文件路径
        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandDi";
        /*vcftools运行文件路径*/
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/";
        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/002_hexaandDi";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/003_tajimaD/001_tajimaD_based100kbWindow_hexa_diploid";
        
        File[] groupFileS = new File(groupDirS).listFiles();
        for (int i = 0; i < groupFileS.length; i++) {
            if (groupFileS[i].isHidden()) {
                groupFileS[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        groupFileS = new File(groupDirS).listFiles();

        //new File(shScriptDirS).mkdir();
        ArrayList<String> perlList = new ArrayList();

        List<String> scriptList = new ArrayList(); //建立一个集合
        try {
            for (int i = 0; i < groupFileS.length; i++) { //对分组进行循环
                String scriptS = new File(ScriptDirS, groupFileS[i].getName().replaceFirst(".txt", "") + "_tajimaD_based100kbWindow.sh").getAbsolutePath();
                scriptList.add(scriptS); //将脚本添加到集合中
                BufferedWriter bw = IOUtils.getTextWriter(scriptS); //开始向脚本中写东西

                String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    if (!chrArr[j].contains("D")) { //说明是不属于D的
//                        String outfileS = new File(outfileDirS, groupFileS[i].getName().replaceFirst(".txt", "_chr" + chrArr[j] + "_based100kbWindow")).getAbsolutePath(); //结果输出的路径
//                        //vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr1D_vmap2.1.vcf --keep /data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/hexaandDi/Ae.tauschii.txt --window-pi 100000 --window-pi-step 50000 --out /data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/001_pi_based100kbWindow_hexa_diploid/chr1D_vmap2.1_Ae.tauschii_ & 
//                        StringBuilder sb = new StringBuilder();
//                        sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath());
//                        sb.append(" --keep ").append(new File(groupFileDirS, groupFileS[i].getName()).getAbsolutePath());
//                        sb.append(" --TajimaD 100000");
//                        sb.append(" --out ").append(outfileS).append("");
//                        bw.write(sb.toString());
//                        bw.newLine();
                    } else {
                        String outfileS = new File(outfileDirS, groupFileS[i].getName().replaceFirst(".txt", "_chr" + chrArr[j] + "_based100kbWindow")).getAbsolutePath(); //结果输出的路径
                        //vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr1D_vmap2.1.vcf --keep /data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/hexaandDi/Ae.tauschii.txt --window-pi 100000 --window-pi-step 50000 --out /data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/001_pi_based100kbWindow_hexa_diploid/chr1D_vmap2.1_Ae.tauschii_ & 
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath());
                        sb.append(" --keep ").append(new File(groupFileDirS, groupFileS[i].getName()).getAbsolutePath());
                        sb.append(" --TajimaD 100000");
                        sb.append(" --out ").append(outfileS).append("");
                        bw.write(sb.toString());
                        bw.newLine();
                    }

                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void mkPiCMDbasedWindow() {
        int windowsize = 100000;
        int windowstep = 50000;

//        //脚本文件路径
//        String ScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/001_script_based100kbWindow/hexaTetra";
//        //本地分组文件路径
//        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandTetra/";
//        /*vcftools运行文件路径*/
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/";
//        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/001_hexaandTetra";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/002_pi_based100kbWindow_hexa_tetraploid";

//        //脚本文件路径
//        String ScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/001_script_based100kbWindow/hexaDi";
//        //本地分组文件路径
//        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandDi";
//        /*vcftools运行文件路径*/
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/";
//        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/002_hexaandDi";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/001_pi_based100kbWindow_hexa_diploid";

//        //脚本文件路径
//        String ScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/003_script_based100kbWindow_byPloidy/hexaandTetra";
//        //本地分组文件路径
//        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/007_byPloidy_forPop/hexaandTetra";
//        /*vcftools运行文件路径*/
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/";
//        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/003_byPloidy_forPop/hexaandTetra";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/003_pi_based100kbWindow_byPloidy/001";
        
        //脚本文件路径
        String ScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/002_Pi/003_script_based100kbWindow_byPloidy/hexaandDi";
        //本地分组文件路径
        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/007_byPloidy_forPop/hexaandDi";
        /*vcftools运行文件路径*/
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/";
        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/003_byPloidy_forPop/hexaandDi";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/003_pi_based100kbWindow_byPloidy/001";


        File[] groupFileS = new File(groupDirS).listFiles();
        for (int i = 0; i < groupFileS.length; i++) {
            if (groupFileS[i].isHidden()) {
                groupFileS[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        groupFileS = new File(groupDirS).listFiles();

        //new File(shScriptDirS).mkdir();
        ArrayList<String> perlList = new ArrayList();

        List<String> scriptList = new ArrayList(); //建立一个集合
        try {
            for (int i = 0; i < groupFileS.length; i++) { //对分组进行循环
                String scriptS = new File(ScriptDirS, groupFileS[i].getName().replaceFirst(".txt", "") + "_Pi_based100kbWindow.sh").getAbsolutePath();
                scriptList.add(scriptS); //将脚本添加到集合中
                BufferedWriter bw = IOUtils.getTextWriter(scriptS); //开始向脚本中写东西

                String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    if (!chrArr[j].contains("D")) { //说明是不属于D的
                        String outfileS = new File(outfileDirS, groupFileS[i].getName().replaceFirst(".txt", "_chr" + chrArr[j] + "_based100kbWindow")).getAbsolutePath(); //结果输出的路径
//                        //vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr1D_vmap2.1.vcf --keep /data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/hexaandDi/Ae.tauschii.txt --window-pi 100000 --window-pi-step 50000 --out /data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/001_pi_based100kbWindow_hexa_diploid/chr1D_vmap2.1_Ae.tauschii_ & 
//                        StringBuilder sb = new StringBuilder();
//                        sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath());
//                        sb.append(" --keep ").append(new File(groupFileDirS, groupFileS[i].getName()).getAbsolutePath());
//                        sb.append(" --window-pi 100000").append(" --window-pi-step 50000 ");
//                        sb.append(" --out ").append(outfileS).append("");
//                        bw.write(sb.toString());
//                        bw.newLine();
                    } else {
                        String outfileS = new File(outfileDirS, groupFileS[i].getName().replaceFirst(".txt", "_chr" + chrArr[j] + "_based100kbWindow")).getAbsolutePath(); //结果输出的路径
                        //vcftools --vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef/chr1D_vmap2.1.vcf --keep /data4/home/aoyue/vmap2/analysis/021_popGen/000_group/002_forPi/hexaandDi/Ae.tauschii.txt --window-pi 100000 --window-pi-step 50000 --out /data4/home/aoyue/vmap2/analysis/021_popGen/002_Pi/001_pi_based100kbWindow_hexa_diploid/chr1D_vmap2.1_Ae.tauschii_ & 
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath());
                        sb.append(" --keep ").append(new File(groupFileDirS, groupFileS[i].getName()).getAbsolutePath());
                        sb.append(" --window-pi 100000").append(" --window-pi-step 50000 ");
                        sb.append(" --out ").append(outfileS).append("");
                        bw.write(sb.toString());
                        bw.newLine();
                    }

                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 目的：将所有sh文本的cmd合并成一个文件。
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergeSh(String infileDirS, String outfileS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sh");
        Arrays.sort(fs);

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = IOUtils.getTextReader(infileS);
                String temp = null; //not read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 目的：将所有sh文本的cmd合并成一个文件,并进行拆分。
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergeandsplitSh(String infileDirS, String outfileS,int numfile,int numcmd) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sh");
        Arrays.sort(fs);

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            //读正文部分
            int cnttotal = 0;
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = IOUtils.getTextReader(infileS);
                String temp = null; //not read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            System.out.println("Total cmd is" + "\t" + cnttotal);
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        String parentS = new File(outfileS).getParent();
        new File(parentS, "splitScript").mkdirs();
        String outfileDirS = new File(parentS, "splitScript").getAbsolutePath();
        String shfileS = new File(parentS, "sh_split.sh").getAbsolutePath();
        String nameprefix = new File(outfileS).getName().replaceFirst(".sh", "_");

        try {
            String[] outS = new String[numfile];
            BufferedReader br = IOUtils.getTextReader(outfileS);
            BufferedWriter[] bw = new BufferedWriter[numfile];
            for (int i = 0; i < outS.length; i++) {
                String num = PStringUtils.getNDigitNumber(3, i + 1);
                outS[i] = new File(outfileDirS, nameprefix + num + ".sh").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outS[i]);
                String temp;
                for (int j = 0; j < numcmd; j++) {
                    if ((temp = br.readLine()) != null) {
                        bw[i].write(temp);
                        bw[i].newLine();
                    }
                }
                bw[i].flush();
                bw[i].close();
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            File[] fss = new File(outfileDirS).listFiles();
            fss = IOUtils.listFilesEndsWith(fss, ".sh");
            Arrays.sort(fss);
            BufferedWriter bw = IOUtils.getTextWriter(shfileS);
            for (int i = 0; i < fss.length; i++) {
                bw.write("sh " + fss[i].getName() + " > log_" + fss[i].getName().split(".sh")[0] + ".txt 2>&1 &");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkFstCommandbasedwinndow() {
//        //local path
//        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandTetra";
//        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandDi";
//        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/002_script_based10kbWindow/hexaTetra";
//        String shScript2DirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/002_script_based10kbWindow/hexaDi";
//
//        //HPC Path
//        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/000_hexaploidandTetraploid";
//        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/001_hexaploidandDiploid";
//
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
//        String outputDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/004_fst_based100kbWindow_hexa_tetraploid";
//        String output2DirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/005_fst_based100kbWindow_hexa_diploid";
        
        
        //local path
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/007_byPloidy_forPop/hexaandTetra_forPi";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/007_byPloidy_forPop/hexaandDi_forPi";
        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/004_script_based10KbWindow_byPloidy/hexaTetra";
        String shScript2DirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/004_script_based10KbWindow_byPloidy/hexaDi";

        //HPC Path
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/003_byPloidy_forPop/hexaandTetra_forPi";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/003_byPloidy_forPop/hexaandDi_forPi";

        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/006_fst_based100kbWindow_hexaandTetra_byPloidy";
        String output2DirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/007_fst_based100kbWindow_hexaandDi_byPloidy";
        
        
        

        File[] group1FileS = new File(groupHexaandTetraDirS).listFiles();
        for (int i = 0; i < group1FileS.length; i++) {
            if (group1FileS[i].isHidden()) {
                group1FileS[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        group1FileS = new File(groupHexaandTetraDirS).listFiles();

        File[] group2FileS = new File(groupHexaandDiDirS).listFiles();
        for (int i = 0; i < group2FileS.length; i++) {
            if (group2FileS[i].isHidden()) {
                group2FileS[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        group2FileS = new File(groupHexaandDiDirS).listFiles();



        new File(shScriptDirS).mkdirs();
        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。
        for (int i = 0; i < group1FileS.length - 1; i++) {
            String group1 = group1FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group1FileS.length; j++) {
                String group2 = group1FileS[j].getName().replace(".txt", "");
                try {
                    String scriptS = new File(shScriptDirS, group1 + "VS" + group2 + ".sh").getAbsolutePath(); //写入的是命令，每个文件包含n条染色体的命令。
                    BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                    for (int k = 1; k < 8; k++) {
                        String[] chr = {k + "A", k + "B"};
                        for (int l = 0; l < chr.length; l++) {
                            String outfileS = new File(outputDirS, group1 + "VS" + group2 + "_chr" + chr[l] + "_").getAbsolutePath();
                            StringBuilder sb = new StringBuilder();
                            sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chr[l] + "_vmap2.1.vcf").getAbsolutePath());
                            sb.append(" --weir-fst-pop ").append(new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath()).append(" --weir-fst-pop ").append(new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath());
                            sb.append(" --fst-window-size 100000").append(" --fst-window-step 50000 ");
                            sb.append(" --out ").append(outfileS).append("");
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }

                    bw.flush();
                    bw.close();
                    perlList.add(new File(scriptS).getName());

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String group1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String group2 = group2FileS[j].getName().replace(".txt", "");
                try {
                    String scriptS = new File(shScript2DirS, group1 + "VS" + group2 + ".sh").getAbsolutePath(); //写入的是命令，每个文件包含n条染色体的命令。
                    BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                    //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                    for (int k = 1; k < 8; k++) {
                        String[] chr = {k + "D"};
                        for (int l = 0; l < chr.length; l++) {
                            String outfileS = new File(output2DirS, group1 + "VS" + group2 + "_chr" + chr[l] + "_").getAbsolutePath();
                            StringBuilder sb = new StringBuilder();
                            sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chr[l] + "_vmap2.1.vcf").getAbsolutePath());
                            sb.append(" --weir-fst-pop ").append(new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath()).append(" --weir-fst-pop ").append(new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath());
                            sb.append(" --fst-window-size 100000").append(" --fst-window-step 50000 ");
                            sb.append(" --out ").append(outfileS).append("");
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }

                    bw.flush();
                    bw.close();
                    perlList.add(new File(scriptS).getName());

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

    }

    /**
     * 求每个值的平均值，做成表格
     */
    public void fstTable() {
//        String groupDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/000_hexaploidandTetraploid";
//        String fstDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/001_fst_basedSNP";
//        String outfileS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/003_fstTable/fstTable_chr001.txt";

        String groupDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/001_hexaploidandDiploid";
        String fstDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/002_fst_basedSNP_hexa_diploid";
        String outfileS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/003_fstTable/fstTable_chr005.txt";

        File[] fs = new File(groupDirS).listFiles(); //求分组文件
        String[] groups = new String[fs.length]; //一个分组文件是一个group
        for (int i = 0; i < groups.length; i++) {
            groups[i] = fs[i].getName().replaceFirst(".txt", "");
        }
        Arrays.sort(groups);
        double[][] fstValues = new double[groups.length][groups.length]; //两两分组
        for (int i = 0; i < groups.length - 1; i++) {
            for (int j = i + 1; j < groups.length; j++) {
//                String infileS = groups[i]+"VS"+groups[j]+"_chr001_.weir.fst"; //先找到其名字
                String infileS = groups[i] + "VS" + groups[j] + "_chr005_.weir.fst"; //先找到其名字
                infileS = new File(fstDirS, infileS).getAbsolutePath(); //再new新的文件，包含绝对路径
                if (!(new File(infileS).exists())) { //判断该文件是否存在
//                    infileS = groups[j]+"VS"+groups[i]+"_chr001_.weir.fst";
                    infileS = groups[j] + "VS" + groups[i] + "_chr005_.weir.fst";
                    infileS = new File(fstDirS, infileS).getAbsolutePath();
                }
                try {
                    BufferedReader br = IOUtils.getTextReader(infileS); //未压缩
                    String temp = br.readLine(); //读表头
                    TDoubleArrayList vList = new TDoubleArrayList();
                    while ((temp = br.readLine()) != null) {
                        List<String> l = PStringUtils.fastSplit(temp);
                        String t = l.get(2);
                        if (t.startsWith("-n") || t.startsWith("n")) {
                            continue;
                        }
                        double v = Double.valueOf(t);
                        if (v < 0) {
                            v = 0;
                        }
                        vList.add(v);
                    }
                    br.close();
                    double[] v = vList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(v);
                    fstValues[i][j] = d.getMean();
                    System.out.println(infileS + String.format("%.3f", fstValues[i][j]));
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Pupoluation");
            for (int i = 0; i < groups.length; i++) { //开始打表头
                sb.append("\t").append(groups[i]);
            }
            bw.write(sb.toString());
            bw.newLine();

            for (int i = 0; i < groups.length; i++) {
                sb = new StringBuilder(groups[i]); //第一列都是分组信息
                for (int j = 0; j < groups.length; j++) {
                    if (fstValues[i][j] == 0) { //说明2个小组没有分化，或者没有比较
                        if (i == j) { //两者相等
                            sb.append("\tNA");
                        } else { //两个group不相等，但是fst值等于0
                            sb.append("\t"); //跳过，说明2者没有进行比较
                        }
                    } else { //两个group不相等，fst值不等于0
                        sb.append("\t").append(fstValues[i][j]);
                    }
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkFstCommandbasedSNP() {
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/000_groups/003_forFst/hexaandDi";
        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/001_Fst/000_scriptsSNPbased";

        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/000_hexaploidandTetraploid";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/000_group/001_hexaploidandDiploid";

        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII";
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/001_fst_basedSNP";
        String output2DirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/001_Fst/002_fst_basedSNP_hexa_diploid";

        File[] group1FileS = new File(groupHexaandTetraDirS).listFiles();
        for (int i = 0; i < group1FileS.length; i++) {
            if (group1FileS[i].isHidden()) {
                group1FileS[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        group1FileS = new File(groupHexaandTetraDirS).listFiles();

        File[] group2FileS = new File(groupHexaandDiDirS).listFiles();
        for (int i = 0; i < group2FileS.length; i++) {
            if (group2FileS[i].isHidden()) {
                group2FileS[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        group2FileS = new File(groupHexaandDiDirS).listFiles();

        new File(shScriptDirS).mkdir();
        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。
//        for (int i = 0; i < group1FileS.length - 1; i++) {
//            String group1 = group1FileS[i].getName().replace(".txt", "");
//            for (int j = i + 1; j < group1FileS.length; j++) {
//                String group2 = group1FileS[j].getName().replace(".txt", "");
//                try {
//                    String scriptS = new File(shScriptDirS, group1 + "VS" + group2 + ".sh").getAbsolutePath(); //写入的是命令，每个文件包含n条染色体的命令。
//                    BufferedWriter bw = IOUtils.getTextWriter(scriptS);
//                    //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
//                    int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
//                    Arrays.sort(db);
//                    for (int k = 1; k < 43; k++) {
//                        String chr = PStringUtils.getNDigitNumber(3, k);
//                        String outfileS = new File (outputDirS, group1 + "VS" + group2 + "_chr" + chr + "_").getAbsolutePath();
//                        if (Arrays.binarySearch(db, k) < 0) { //是属于AB的
//                            StringBuilder sb = new StringBuilder();
//                            sb.append("vcftools --vcf ").append(new File(infileDirS,"chr" + chr + "_vmap2.1.vcf").getAbsolutePath()).append(" --weir-fst-pop ").append(new File(group1FileDirS,group1FileS[i].getName()).getAbsolutePath()).append(" --weir-fst-pop ").
//                                    append(new File(group1FileDirS,group1FileS[j].getName()).getAbsolutePath()).append(" --out ").append(outfileS).append("");
//                            bw.write(sb.toString());
//                            bw.newLine();
//                        }
//                    }
//                    bw.flush();
//                    bw.close();
//                    perlList.add(new File(scriptS).getName());
//
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//            }
//        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String group1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String group2 = group2FileS[j].getName().replace(".txt", "");
                try {
                    String scriptS = new File(shScriptDirS, group1 + "VS" + group2 + ".sh").getAbsolutePath(); //写入的是命令，每个文件包含n条染色体的命令。
                    BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                    //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                    int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
                    Arrays.sort(db);
                    for (int k = 1; k < 43; k++) {
                        String chr = PStringUtils.getNDigitNumber(3, k);
                        String outfileS = new File(output2DirS, group1 + "VS" + group2 + "_chr" + chr + "_").getAbsolutePath();
                        if (Arrays.binarySearch(db, k) > -1) { //是属于D的
                            StringBuilder sb = new StringBuilder();
                            sb.append("vcftools --vcf ").append(new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath()).append(" --weir-fst-pop ").append(new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath()).append(" --weir-fst-pop ").
                                    append(new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath()).append(" --out ").append(outfileS).append("");
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                    bw.flush();
                    bw.close();
                    perlList.add(new File(scriptS).getName());

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        /**
         * ****** 将15个组合的脚本命令统一执行，有一个统一执行的脚本 sh_fst.sh ******
         */
        try {
            String outfileS = new File(shScriptDirS, "sh_fst.sh").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < perlList.size(); i++) {
                StringBuilder sb = new StringBuilder("sh ");
                sb.append(perlList.get(i))
                        .append("");;
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void getTaxaSet() {
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/abd/002_taxaRefBam.ABDgenome.manual.addNAFU.txt";
            String outfileS = "/Users/Aoyue/Documents/BreadWheat.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            Set<String> s = new HashSet<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                l = PStringUtils.fastSplit(temp);
                s.add(l.get(0));
                System.out.println(l.get(0));
            }
            br.close();
            System.out.println(s.size() + " taxa");

            String[] sS = s.toArray(new String[s.size()]);
            Arrays.sort(sS);
            for (int i = 0; i < s.size(); i++) {
                bw.write(sS[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
