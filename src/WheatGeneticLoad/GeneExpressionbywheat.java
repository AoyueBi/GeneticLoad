package WheatGeneticLoad;

import AoUtils.*;
import AoUtils.Triads.Triadsgenes;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.lang.ArrayUtils;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class GeneExpressionbywheat {
    public GeneExpressionbywheat(){
//        this.addTPM();
//        this.getCSgeneSummary();
//        this.getCSgeneSummary_usingTriadsClass();

//        this.geneExpressionbyDevelopment();
//        this.geneExpressionbyRoot_fromJun();
//        this.geneExpressionbyColeoptile_fromXiaohan();

        /**
         * tissue express breadth
         */
        this.getTissueBreadth();
        this.getwindowDistrbution();

        this.mergegeneExpression_fromJun();
//        this.calBreadth_onRootandcoleotiple();
//        this.getwindowDistrbution_general();

    }


    /**
     * 将每个基因对应的 value 进行 window 扫描 500 个值 100/500=0.2 bin 宽度
     */
    public void getwindowDistrbution_general(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_Hexaploid_root_geneExpression.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/003_window/001_Hexaploid_root_geneExpression.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/003_Hexaploid_coleoptile_geneExpression.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/003_window/001_Hexaploid_coleoptile_geneExpression.txt";

        AoFile.readheader(infileS);
        String[] chrS = AoFile.getStringArraybySet(infileS,0);
        TDoubleArrayList[] posList = new TDoubleArrayList[chrS.length];
        TDoubleArrayList[] valueList = new TDoubleArrayList[chrS.length];
        for (int i = 0; i < chrS.length; i++) { //记得初始化
            posList[i] = new TDoubleArrayList();
            valueList[i] = new TDoubleArrayList();
        }
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0);
                double posScale = Double.parseDouble(l.get(5)); //**************** 需要修改
                double value = Double.parseDouble(l.get(3)); //**************** 需要修改
                int index = Arrays.binarySearch(chrS,chr);
                if (index <0){
                    System.out.println(chr + "\t" + posScale);
                }
                posList[index].add(posScale);
                valueList[index].add(value);
            }
            br.close();

            List[][] output = new List[chrS.length][];
            for (int i = 0; i < chrS.length; i++) {
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,1,0.5);
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,2,1);
                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,4,2);

            }

            bw.write("Chr\tPos_scale\tCount\tAveExpression");
            bw.newLine();
            for (int i = 0; i < chrS.length; i++) { //第一层循环是1A到7D
                List<String>[] out = output[i];
                for (int j = 0; j < out[0].size(); j++) { //第二层循环是
                    String chr = chrS[i];
                    String pos = out[0].get(j);
                    String count = out[1].get(j);
                    String value = out[2].get(j);
                    bw.write(chr + "\t" + pos + "\t" + count + "\t" + value);
                    bw.newLine();
                }
            }
            int a=3;
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * 将每个基因对应的 value 进行 window 扫描 500 个值 100/500=0.2 bin 宽度
     */
    public void getwindowDistrbution_onRootColeotiple(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/007_Hexaploid_root_coleotiple_geneExpression_nofilter_ExpressionBreadth.txt";
        AoFile.readheader(infileS);
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/008_Hexaploid_root_coleotiple_geneExpression_nofilter_ExpressionBreadth_window.txt";
        String[] chrS = AoFile.getStringArraybySet(infileS,0);
        TDoubleArrayList[] posList = new TDoubleArrayList[chrS.length];
        TDoubleArrayList[] valueList = new TDoubleArrayList[chrS.length];
        for (int i = 0; i < chrS.length; i++) { //记得初始化
            posList[i] = new TDoubleArrayList();
            valueList[i] = new TDoubleArrayList();
        }
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0);
                double posScale = Double.parseDouble(l.get(2));
                double value = Double.parseDouble(l.get(8));
                int index = Arrays.binarySearch(chrS,chr);
                if (index <0){
                    System.out.println(chr + "\t" + posScale);
                }
                posList[index].add(posScale);
                valueList[index].add(value);
            }
            br.close();

            List[][] output = new List[chrS.length][];
            for (int i = 0; i < chrS.length; i++) {
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,1,0.5);
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,2,1);
                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,4,2);

            }

            bw.write("Chr\tPos_scale\tCount\tExpressionBreadth");
            bw.newLine();
            for (int i = 0; i < chrS.length; i++) { //第一层循环是1A到7D
                List<String>[] out = output[i];
                for (int j = 0; j < out[0].size(); j++) { //第二层循环是
                    String chr = chrS[i];
                    String pos = out[0].get(j);
                    String count = out[1].get(j);
                    String value = out[2].get(j);
                    bw.write(chr + "\t" + pos + "\t" + count + "\t" + value);
                    bw.newLine();
                }
            }
            int a=3;
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    private void calBreadth_onRootandcoleotiple(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/006_Hexaploid_root_coleotiple_geneExpression_nofilter.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/007_Hexaploid_root_coleotiple_geneExpression_nofilter_ExpressionBreadth.txt";
        try {

            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header + "\tExpressionBreadth");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cntBreadth = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                double root = Double.parseDouble(l.get(4));
                double coleotiple = Double.parseDouble(l.get(6));
                if (root < 0.5 && coleotiple < 0.5) continue;
                if (root < 0.5 && coleotiple >= 0.5){
                    cntBreadth = 1;
                }
                if (root >= 0.5 && coleotiple < 0.5){
                    cntBreadth = 1;
                }
                if (root >= 0.5 && coleotiple >= 0.5){
                    cntBreadth = 2;
                }
                bw.write(temp + "\t" + cntBreadth);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * step 1: 将根和胚芽鞘的数据合并起来，进行组织表达宽度分析，不过滤，都是53169个基因
     *
     * Chr\tPos_start\tGene\tAveExpression\tSD
     */
    public void mergegeneExpression_fromJun(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/source/7_root_nor_countResult.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/004_Hexaploid_root_geneExpression_nofilter.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/source/3_nor_countResult.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/005_Hexaploid_cpleoptile_geneExpression_nofilter.txt";

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        gf.sortGeneByName();
        Triadsgenes tg = new Triadsgenes();
        System.out.println(tg.getGeneNum());
        AoFile.readheader(infileS);


        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Chr\tPos_start\tGene\tAveExpression\tSD\tPos_scaled");
            bw.newLine();
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt =0;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String gene = l.get(0);
                boolean bl = tg.ifTriads(gene);
                if (!bl)continue; // ******** filter gene which is not in triad
                boolean blexpression = tg.ifExpressedBasedGene(gene);
                if (!blexpression)continue; // ******** filter gene which is not expressed
                cnt++;
                TDoubleArrayList tpmList = new TDoubleArrayList();
                for (int i = 1; i < l.size(); i++) {
                    tpmList.add(Double.parseDouble(l.get(i)));
                }
                String ave = AoMath.getRelativeMean(tpmList);
                String sd = AoMath.getStandardDeviation(tpmList);
//                if (Double.parseDouble(ave) < 0.5)continue; //******** filter TPM ave less than 0.5
                //get gene chr pos start
                int index = gf.getGeneIndex(gene);
                int chr = gf.getGeneChromosome(index);
                int pos = gf.getGeneStart(index);
                String chromosome = RefV1Utils.getChromosome(chr,pos);
                int posOnchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                String posScaled = WheatUtils.getScaledPos(chromosome,posOnchromosome);
                bw.write(chromosome + "\t" + posOnchromosome + "\t" + gene + "\t" + ave + "\t" + sd + "\t" + posScaled);
                bw.newLine();
            }
            System.out.println(cnt + " genes kept");
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }



    /**
     * 将每个基因对应的 value 进行 window 扫描 500 个值 100/500=0.2 bin 宽度
     */
    public void getwindowDistrbution(){

//        String infileS = "";
//        String outfileS = "";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_Azhurnaya_tissue_breadth_miniVersion.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/003_Azhurnaya_tissue_breadth_window.txt";


//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/005_Azhurnaya_tissue_breadth_100Kgene_miniVersion.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/006_Azhurnaya_tissue_breadth_100Kgene_window.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/002_Azhurnaya_tissue_breadth_100Kgenet_1TPM_miniVersion.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/003_Azhurnaya_tissue_breadth_100Kgenet_1TPM_miniVersion_window.txt";

//                String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/002_Azhurnaya_tissue_breadth_100Kgenet_50TPM_miniVersion.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/002_Azhurnaya_tissue_breadth_100Kgenet_50TPM_miniVersion_window.txt";

//                String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/002_Azhurnaya_tissue_breadth_100Kgenet_100TPM_miniVersion.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/002_Azhurnaya_tissue_breadth_100Kgenet_100TPM_miniVersion_window.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/005_Azhurnaya_tissue_breadth_100Kgene_miniVersion.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/007_Azhurnaya_tissue_breadth_100Kgene_miniVersion_windowBypos.txt";

        AoFile.readheader(infileS);

        String[] chrS = AoFile.getStringArraybySet(infileS,0);
        TDoubleArrayList[] posList = new TDoubleArrayList[chrS.length];
        TDoubleArrayList[] valueList = new TDoubleArrayList[chrS.length];
        for (int i = 0; i < chrS.length; i++) { //记得初始化
            posList[i] = new TDoubleArrayList();
            valueList[i] = new TDoubleArrayList();
        }
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0);
                double posScale = Double.parseDouble(l.get(1));
//                double posScale = Double.parseDouble(l.get(2));
                double value = Double.parseDouble(l.get(4));
                int index = Arrays.binarySearch(chrS,chr);
                if (index <0){
                    System.out.println(chr + "\t" + posScale);
                }
                posList[index].add(posScale);
                valueList[index].add(value);
            }
            br.close();

            List[][] output = new List[chrS.length][];
            for (int i = 0; i < chrS.length; i++) {
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,1,0.5);
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,2,1);
//                output[i]= Bin.windowstep_posAve(posList[i],valueList[i],100,4,2);
                output[i]= Bin.windowstep_posAve2(posList[i],valueList[i],4,2);


            }

//            bw.write("Chr\tPos_scale\tCount\tExpressionBreadth");
            bw.write("Chr\tPos\tCount\tExpressionBreadth");
            bw.newLine();
            for (int i = 0; i < chrS.length; i++) { //第一层循环是1A到7D
                List<String>[] out = output[i];
                for (int j = 0; j < out[0].size(); j++) { //第二层循环是
                    String chr = chrS[i];
                    String pos = out[0].get(j);
                    String count = out[1].get(j);
                    String value = out[2].get(j);
                    bw.write(chr + "\t" + pos + "\t" + count + "\t" + value);
                    bw.newLine();
                }
            }
            int a=3;
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据azhurnaya的22个组织，求每个组织的平均表达量和方差，小于0.5就算这个组织不表达。
     * 2次测试：第一是只看1：1：1的基因；第二是看pgf中的107891个基因；
     * 表格形式：
     * Chr\tPos_start\tGene\tT1_AveExpression(SD)\tT2_AveExpression\tT3_AveExpression.... \tT22_AveExpression（ave）\tTissueBreadth
     */
    public void getTissueBreadth(){

        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/wheatRNAScience/download/Development_tpm.tsv.gz";
        String AzhurnayaInfoS = "/Users/Aoyue/Documents/Data/wheat/article/wheatRNAScience/myself/Azhurnaya_developmentInfo.txt";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/001_Azhurnaya_tissue_breadth.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/004_Azhurnaya_tissue_breadth_100Kgene.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/001_Azhurnaya_tissue_breadth_100Kgenet_1TPM.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/001_Azhurnaya_tissue_breadth_100Kgenet_50TPM.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_tissueBreadth/002_1TPMexpress/001_Azhurnaya_tissue_breadth_100Kgenet_100TPM.txt";

        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        gf.sortGeneByName();
        Triadsgenes tg = new Triadsgenes();
        AoFile.readheader(infileS);

        /**
         * get tissue's sample
         */
        String[] tissueArray = AoFile.getStringArraybySet(AzhurnayaInfoS,6);
        for (int i = 0; i < tissueArray.length; i++) {
            System.out.println(tissueArray[i]);
        }
        Arrays.sort(tissueArray);
        /**
         * sample <-> tissue HashMap
         */
        HashMap<String,String> hmSampleTissue = AoFile.getHashMapStringKey(AzhurnayaInfoS,0,6);
        /**
         * each tissue contains several samples, we need to know the sample index
         */
        List<Integer>[] tissueIndexList = new List[tissueArray.length];
        for (int i = 0; i < tissueArray.length; i++) {
            tissueIndexList[i] = new ArrayList<>();
        }

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            //deal with header
            String header = br.readLine();
            l = PStringUtils.fastSplit(header);
            for (int i = 1; i < l.size(); i++) { //从 index 1列开始
                String sampleID = l.get(i);
                String tissue = hmSampleTissue.get(sampleID);
                int index = Arrays.binarySearch(tissueArray,tissue);
                if (index < 0){
                    System.out.println(sampleID + "\tin the expression file is not in the Tissue arrays");
                    continue;
                }
                tissueIndexList[index].add(i); //获取索引
            }

            bw.write("Chr\tPos\tPos_scale\tGene");
            for (int i = 0; i < tissueArray.length; i++) {
                bw.write("\t" + tissueArray[i]+ "_AveExpression" + "\t" + tissueArray[i]+"_SD");
            }
            bw.write("\t" + "ExpressionBreadth");
            bw.newLine();

            // deal with expression values
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String gene = AoString.getv11geneName(l.get(0));

                /**
                 * ############## need to modify identify which genes will you analysis: 1:1:1
                 */
//                boolean bl = tg.ifTriads(gene);
//                if (!bl)continue; // ******** filter gene which is not in triad
//                boolean blexpression = tg.ifExpressedBasedGene(gene);
//                if (!blexpression)continue; // ******** filter gene which is not expressed

                /**
                 * ############## need to modify identify which genes will you analysis: 107891
                 */

                int ind = gf.getGeneIndex(gene);
                if (ind < 0) continue;
                cnt++;

                // 目的： 找到每个组织的样品对应的基因表达量，并加入list中，求平均值方差之类的
                // 已获取： 每个组织对应样本的索引
                // 如何操作？ 第一：进行组织的循环，在每个组织中，进行index的循环，找到对应的index，然后添加到 tissueValueList 中
                TDoubleArrayList[] tissueValueList = new TDoubleArrayList[tissueArray.length];
                for (int i = 0; i < tissueArray.length; i++) {
                    tissueValueList[i] = new TDoubleArrayList();
                }
                String[] ave = new String[tissueArray.length];
                String[] sd = new String[tissueArray.length];
                int breadth = 0;
                for (int i = 0; i < tissueArray.length; i++) { //22 tissues
                    for (int j = 0; j < tissueIndexList[i].size(); j++) { // 第 i 个tissue 中有
                        String value = l.get(tissueIndexList[i].get(j));
                        double valued = Double.parseDouble(value);
                        tissueValueList[i].add(valued);
                    }
                    ave[i] = AoMath.getRelativeMean(tissueValueList[i]);
                    sd[i] = AoMath.getStandardDeviation(tissueValueList[i]);

                } //完成了一个基因所有组织的样本对应的表达量的添加，已经这个组织的平均值和标准差

//                for (int i = 0; i < ave.length; i++) { //判断组织中平均表达量大于0.5的，就算是表达了
//                    if (Double.parseDouble(ave[i]) > 0.5){
//                        breadth++;
//                    }
//                }

                for (int i = 0; i < ave.length; i++) { //判断组织中平均表达量大于 1 的，就算是表达了
                    if (Double.parseDouble(ave[i]) > 100){
                        breadth++;
                    }
                }
                if (breadth == 0) continue; //过滤没有在任何组织表达的基因

                //get gene chr pos start
                int index = gf.getGeneIndex(gene);
                int chr = gf.getGeneChromosome(index);
                int pos = gf.getGeneStart(index);
                String chromosome = RefV1Utils.getChromosome(chr,pos);
                int posOnchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                String posScaled = WheatUtils.getScaledPos(chromosome,posOnchromosome);


                bw.write(chromosome + "\t" + posOnchromosome + "\t" + posScaled + "\t" + gene);
                for (int i = 0; i < ave.length; i++) {
                    bw.write("\t" + ave[i] + "\t" + sd[i]);
                }
                bw.write("\t" + breadth);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * step 1: 找所有组织的一个样品的平均基因表达量
     * step 2: 输出文件格式
     * Chr\tPos_start\tGene\tAveExpression\tSD
     */
    public void geneExpressionbyColeoptile_fromXiaohan(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/source/countResult_3_coleoptile_TPM.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/003_Hexaploid_coleoptile_geneExpression.txt";

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        gf.sortGeneByName();
        Triadsgenes tg = new Triadsgenes();
        System.out.println(tg.getGeneNum());
        AoFile.readheader(infileS);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Chr\tPos_start\tGene\tAveExpression\tSD\tPos_scaled");
            bw.newLine();
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt =0;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String gene = l.get(3);
                boolean bl = tg.ifTriads(gene);
                if (!bl)continue; // ******** filter gene which is not in triad
                boolean blexpression = tg.ifExpressedBasedGene(gene);
                if (!blexpression)continue; // ******** filter gene which is not expressed
                cnt++;
                TDoubleArrayList tpmList = new TDoubleArrayList();
                for (int i = 4; i < l.size(); i++) {
                    tpmList.add(Double.parseDouble(l.get(i)));
                }
                String ave = AoMath.getRelativeMean(tpmList);
                String sd = AoMath.getStandardDeviation(tpmList);
                if (Double.parseDouble(ave) < 0.5)continue; //******** filter TPM ave less than 0.5
                //get gene chr pos start
                int index = gf.getGeneIndex(gene);
                int chr = gf.getGeneChromosome(index);
                int pos = gf.getGeneStart(index);
                String chromosome = RefV1Utils.getChromosome(chr,pos);
                int posOnchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                String posScaled = WheatUtils.getScaledPos(chromosome,posOnchromosome);
                bw.write(chromosome + "\t" + posOnchromosome + "\t" + gene + "\t" + ave + "\t" + sd + "\t" + posScaled);
                bw.newLine();
            }
            System.out.println(cnt + " genes kept"); //53196
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * step 1: 找所有组织的一个样品的平均基因表达量
     * step 2: 输出文件格式
     * Chr\tPos_start\tGene\tAveExpression\tSD
     */
    public void geneExpressionbyRoot_fromJun(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/source/7_root_nor_countResult.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/002_Hexaploid_root_geneExpression.txt";


        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/source/3_nor_countResult.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/003_Hexaploid_cpleoptile_geneExpression.txt";

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        gf.sortGeneByName();
        Triadsgenes tg = new Triadsgenes();
        System.out.println(tg.getGeneNum());
        AoFile.readheader(infileS);


        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Chr\tPos_start\tGene\tAveExpression\tSD\tPos_scaled");
            bw.newLine();
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt =0;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String gene = l.get(0);
                boolean bl = tg.ifTriads(gene);
                if (!bl)continue; // ******** filter gene which is not in triad
                boolean blexpression = tg.ifExpressedBasedGene(gene);
                if (!blexpression)continue; // ******** filter gene which is not expressed
                cnt++;
                TDoubleArrayList tpmList = new TDoubleArrayList();
                for (int i = 1; i < l.size(); i++) {
                    tpmList.add(Double.parseDouble(l.get(i)));
                }
                String ave = AoMath.getRelativeMean(tpmList);
                String sd = AoMath.getStandardDeviation(tpmList);
                if (Double.parseDouble(ave) < 0.5)continue; //******** filter TPM ave less than 0.5
                //get gene chr pos start
                int index = gf.getGeneIndex(gene);
                int chr = gf.getGeneChromosome(index);
                int pos = gf.getGeneStart(index);
                String chromosome = RefV1Utils.getChromosome(chr,pos);
                int posOnchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                String posScaled = WheatUtils.getScaledPos(chromosome,posOnchromosome);
                bw.write(chromosome + "\t" + posOnchromosome + "\t" + gene + "\t" + ave + "\t" + sd + "\t" + posScaled);
                bw.newLine();
            }
            System.out.println(cnt + " genes kept");
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * step 1: 找所有组织的一个样品的平均基因表达量
     * step 2: 输出文件格式
     * Chr\tPos_start\tGene\tAveExpression\tSD
     */
    public void geneExpressionbyDevelopment(){

        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/wheatRNAScience/download/Development_tpm.tsv.gz";
        String AzhurnayaInfoS = "/Users/Aoyue/Documents/Data/wheat/article/wheatRNAScience/myself/Azhurnaya_developmentInfo.txt";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        gf.sortGeneByName();
        Triadsgenes tg = new Triadsgenes();
        System.out.println(tg.getGeneNum());
        AoFile.readheader(infileS);


        try {
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/002_testGeneExpressionDistribution/001_Azhurnaya_development_geneExpression.txt";
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Chr\tPos_start\tGene\tAveExpression\tSD\tPos_scaled");
            bw.newLine();
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt =0;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String gene = AoString.getv11geneName(l.get(0));
                boolean bl = tg.ifTriads(gene);
                if (!bl)continue; // ******** filter gene which is not in triad
                boolean blexpression = tg.ifExpressedBasedGene(gene);
                if (!blexpression)continue; // ******** filter gene which is not expressed
                cnt++;
                TDoubleArrayList tpmList = new TDoubleArrayList();
                for (int i = 1; i < l.size(); i++) {
                    tpmList.add(Double.parseDouble(l.get(i)));
                }
                String ave = AoMath.getRelativeMean(tpmList);
                String sd = AoMath.getStandardDeviation(tpmList);
                if (Double.parseDouble(ave) < 0.5)continue; //******** filter TPM ave less than 0.5
                //get gene chr pos start
                int index = gf.getGeneIndex(gene);
                int chr = gf.getGeneChromosome(index);
                int pos = gf.getGeneStart(index);
                String chromosome = RefV1Utils.getChromosome(chr,pos);
                int posOnchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                String posScaled = WheatUtils.getScaledPos(chromosome,posOnchromosome);
                bw.write(chromosome + "\t" + posOnchromosome + "\t" + gene + "\t" + ave + "\t" + sd + "\t" + posScaled);
                bw.newLine();
            }
            System.out.println(cnt + " genes kept");
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }






    /**
     * 选取CS的800个基因，进行相关性分析
     * 去除不表达的基因,并添加 ternary分组和是否共线等信息
     */
    public void getCSgeneSummary_usingTriadsClass(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/001_geneSummary_hexaploid_triad_Removed169_addTPM.txt.gz";
        String csgeneFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/CS_geneSummary_triads_Remove169_wideTable.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/002_geneSummary_CS_triad_testTriadsClass.txt";

        List<String> genesinCS = new ArrayList<>();
        List<String> geneA = AoFile.getStringList(csgeneFileS,5);
        List<String> geneB = AoFile.getStringList(csgeneFileS,6);
        List<String> geneD = AoFile.getStringList(csgeneFileS,7);
        genesinCS.addAll(geneA);genesinCS.addAll(geneB);genesinCS.addAll(geneD);
        HashMap<String,String> hmTriadIDLoadGroup;
        AoFile.readheader(csgeneFileS);
        hmTriadIDLoadGroup = AoFile.getHashMapStringKey(csgeneFileS,0,4);
        Collections.sort(genesinCS);
        Triadsgenes triad = new Triadsgenes();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tTernaryGroup\tIfSyntenic");
            bw.newLine();
            String temp;
            List<String> l;
            int cnt = 0;
            String syntenicState;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(1).split("\\.")[0];
                int index = Collections.binarySearch(genesinCS,gene);
                if (index < 0)continue;
                boolean out = triad.ifExpressedBasedGene(gene); //判断是否表达
                if (!out)continue;
                String triadID = triad.getTriadID(gene);
                String loadGroup = hmTriadIDLoadGroup.get(triadID); //获取 分组信息
                boolean syntenic = triad.ifSyntenicBasedTriadID(triadID);
                if (syntenic) syntenicState = "1";
                else syntenicState = "0";
                bw.write(temp+"\t"+ loadGroup + "\t" + syntenicState);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 选取CS的800个基因，进行相关性分析
     * 去除不表达的基因,并添加 ternary分组和是否共线等信息
     */
    public void getCSgeneSummary(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/001_geneSummary_hexaploid_triad_Removed169_addTPM.txt.gz";
        String csgeneFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/CS_geneSummary_triads_Remove169_wideTable.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/002_geneSummary_CS_triad_Removed169_addTPM_removeNoexpression_removeNOHGDel.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/002_geneSummary_CS_triad_Removed169_addTPM_removeNoexpression.txt.gz";

        List<String> genesinCS = new ArrayList<>();
        HashMap<String,String> hmTriadIDLoadGroup = new HashMap<>();
        HashMap<String,String> hmgeneIfexpressed = new HashMap<>();
        HashMap<String,String> hmgeneIfsyntenic = new HashMap<>();

        AoFile.readheader(csgeneFileS);

        try{
            BufferedReader br = AoFile.readFile(csgeneFileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String triadID = l.get(0);
                String loadGroup = l.get(4);
                String genea = l.get(5);
                String geneb = l.get(6);
                String gened = l.get(7);
                String synteny = l.get(8);
                String expressed = l.get(9);

                /**
                 * 过滤NA 即没有 high del
                 */
                String contribution = l.get(1);
                if (contribution.equals("NA")) continue;

//                if (expressed.equals("FALSE"))continue;
                genesinCS.add(genea);genesinCS.add(geneb);genesinCS.add(gened);
                hmgeneIfexpressed.put(genea,expressed);
                hmgeneIfexpressed.put(geneb,expressed);
                hmgeneIfexpressed.put(gened,expressed);
                hmTriadIDLoadGroup.put(genea,loadGroup);
                hmTriadIDLoadGroup.put(geneb,loadGroup);
                hmTriadIDLoadGroup.put(gened,loadGroup);
                hmgeneIfsyntenic.put(genea,synteny);
                hmgeneIfsyntenic.put(geneb,synteny);
                hmgeneIfsyntenic.put(gened,synteny);
            }
            br.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        Collections.sort(genesinCS);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tTernaryGroup\tIfSyntenic");
            bw.newLine();
            String temp;
            List<String> l;
            int cnt = 0;
            String syntenicState = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(1).split("\\.")[0];
                int index = Collections.binarySearch(genesinCS,gene);
                if (index < 0)continue;
                /**
                 * new method
                 */
                String express = hmgeneIfexpressed.get(gene);
                if (express.equals("FALSE"))continue;
                String loadGroup = hmTriadIDLoadGroup.get(gene);
                syntenicState = hmgeneIfsyntenic.get(gene);
                bw.write(temp+"\t"+ loadGroup + "\t" + syntenicState);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 向基因summary文件中添加 TPM 信息，即每个基因的平均表达量，最大值，最小值，标准误
     *
     * step1:
     */
    public  void addTPM (){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/002_geneSummary_hexaploid_triad_Removed169.txt.gz";
        String geneExpressionS = "/Users/Aoyue/Documents/Data/wheat/gene/gene_expression_dontMove/geneExpression.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/001_geneSummary_hexaploid_triad_Removed169_addTPM.txt.gz";

        AoFile.readheader(infileS);
        //建库HashMap
        HashMap<String,Double> hmgeneMean = AoFile.getHashMapdoubleValue(geneExpressionS,0,5);
        HashMap<String,Double> hmgeneSd = AoFile.getHashMapdoubleValue(geneExpressionS,0,6);
        HashMap<String,Double> hmgeneRsd = AoFile.getHashMapdoubleValue(geneExpressionS,0,7);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tTPM_mean\tTPM_sd\tTPM_rsd");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(1).split("\\.")[0];
                double meanTPM = hmgeneMean.get(gene);
                double sdTPM = hmgeneSd.get(gene);
                double rsdTPM = hmgeneRsd.get(gene);
                bw.write(temp+"\t"+String.format("%.4f",meanTPM) + "\t" + String.format("%.4f",sdTPM)+ "\t" + String.format("%.4f",rsdTPM));
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


}
