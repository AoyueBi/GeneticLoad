package PopulationAnalysis;

import AoUtils.Bin;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class Fst {
    public Fst() {
//        this.mkFstCommandbasedSNP();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/001_scriptSNPbased/fst_basedSNP_20200205.sh",7,23);

//        this.mkFstTable("/Users/Aoyue/Documents/po1_VS_pop2_chr1A.txt","/Users/Aoyue/Documents/out.txt");

        this.scriptMkFstTable();
    }


    public void mkFstCommandbasedwinndow() {

        //local file： one: group two: script
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";
        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_scriptSNPbased";

        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS = "";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        new File(shScriptDirS).mkdirs(); //无论在不在，都可以创建该文件夹

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }
    }


    /**
     * 进行批量计算
     *
     */
    public void scriptMkFstTable() {

        //local file： one: group two: script
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";
        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_scriptSNPbased";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/002_fst_basedSNP";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/003_fst_table";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/log/001";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        new File(shScriptDirS).mkdirs(); //无论在不在，都可以创建该文件夹

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".weir.fst").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    String logS = new File(logDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    System.out.println("java -jar 036_mkFstTable.jar " + infileS + " " + outfileS + " > " + logS + " 2>&1");
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".weir.fst").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    String logS = new File(logDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    System.out.println("java -jar 036_mkFstTable.jar " + infileS + " " + outfileS + " > " + logS + " 2>&1");
                }
            }
        }
    }

    public void scriptMkFstTable2(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/002_fst_basedSNP";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/003_fst_table";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/log/001";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            if (infileS.endsWith(".txt")) {
                outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_RH_2Mwindow_1Mstep.txt.gz").getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_RH_2Mwindow_1Mstep.txt.gz").getAbsolutePath();
            }
            new Bin().calwindowstep_ResidualHeterozygosity(infileS,2000000,1000000,outfileS);
            System.out.println(f.getName() + "\tis completed at " + outfileS);
        });
    }

    /**
     *
     *
     */
    public void mkFstTable(String infileS, String outfileS){
        //pseudo-code
        /**
         * 四倍体：WE DE FT
         * 六倍体：LR CU
         * 二倍体
         * 合计12种组合，两两之间。每种组合分为A B subgenome,每条Asub 分为 1-7条染色体
         *
         * 先跑每一种组合的7条，写成1A 2A 3A 4A 5A 6A 7A
         * 只要知道infileS,就能得出 chr group value , 就能输出到 outfileS 中；
         */

        try {
            String name = new File(infileS).getName();
            String chr = this.getChr(infileS);
            String pop1pop2 = name.substring(0,name.indexOf("_chr"));
            String value = this.getMean(infileS);

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tPOP1_VS_POP2\tMEAN_FST");
            bw.newLine();
            bw.write(chr+ "\t" + pop1pop2 + "\t" + value);
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据文件的名字获取chr信息，不进行里面的读取
     *
     * @param infileS
     * @return
     */
    public String getChr(String infileS){
        String out = null;

        String name = new File(infileS).getName();
        int index = name.indexOf("chr");
        int index1 = index+3;
        int index2 = index+5;
        out = name.substring(index1,index2);
        System.out.println(out);
        return out;
    }


    /**
     * get mean from a file， is the value equals -nan, continue
     *
     * @param infileS
     * @return
     */
    public String getMean(String infileS){
        String out = null;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            TDoubleArrayList vList = new TDoubleArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String t = l.get(2);
                if (t.startsWith("-n") || t.startsWith("n")) { //pos中 fst值为 -nan 的全部去掉
                    continue;
                }
                double v = Double.parseDouble(t);
                vList.add(v);
                cnt++; //有 fst值的位点数目
            }
            br.close();
            double[] v = vList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(v);
            double m = d.getMean();
            out = String.format("%.4f", m);
            System.out.println(infileS + " with mean fst " + out + " is completed");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return out;
    }

    /**
     *
     */
    public void mkFstCommandbasedSNP() {

        //local file： one: group two: script
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";
        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_scriptSNPbased";

        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/002_fst_basedSNP";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        new File(shScriptDirS).mkdirs(); //无论在不在，都可以创建该文件夹

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outputDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outputDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }
    }

}



