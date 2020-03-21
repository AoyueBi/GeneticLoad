/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import pgl.format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.utils.IOUtils;
import pgl.utils.PArrayUtils;
import pgl.utils.PStringUtils;
import pgl.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class Bin {

    public Bin() {
        //this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/002_calMAF","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/003_bintable","25","0.5");

//    this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/009_calMAF_newData","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/010_bintable","25","0.5");
//    this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/012_calMAF_bySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/013_bintable_bySub","25","0.5");
//        this.getDAFtable();
    }


    /**
     * 根据 chr pos 和 value 来确定,返回 每个window内的变异个数，以及 残余杂合度 Residual Heterozygosity
     *
     */
    public void calwindowstep_ResidualHeterozygosity(String infileS, int window, int step, String outfileS){
        // 1.将pos转为为list,找到最大值，根据最大值确定bin的数目
        // 2.建立count数目和list数组，对pos进行循环，找到每个bin的左边的数目和value的集合，求这个集合的平均值，最大值，方差等
        // 3.输出，每个bin的值
        //************************

        //先求最大值,即染色体长度
        RowTable<String> t = new RowTable(infileS);
        int chrlength = Integer.valueOf(t.getCell(t.getRowNumber() - 1, 1)); //t.getRowNumber()是文件的行数，不包括header。 这里getCell得到的是索引         //染色体的长度是最后一行pos的位置
        String chr = t.getCell(0, 0);
        System.out.println(new File(infileS).getName().substring(0,5) + " length is " + chrlength);

        int[][] bound = this.initializeWindowStep(chrlength, window,step);

        int count[] = new int[bound.length]; //查看每个bin里面的变异个数
        List<String>[] value = new ArrayList[bound.length]; //每个Bin里面的值的集合 List

        int[] boundright = new int[bound.length]; //只看左边的bound
        int[] boundleft = new int[bound.length]; //右边的bound
        for (int i = 0; i < bound.length; i++) { //每个bound的左边
            boundleft[i] = bound[i][0];
            boundright[i] = bound[i][1];
            value[i] = new ArrayList(); //对每一个bin中的List进行初始化
        }


        for (int i = 0; i < t.getRowNumber(); i++) {
            int pos = Integer.parseInt(t.getCell(i,1));
            String va = t.getCell(i,2);

            int indexleft = Arrays.binarySearch(boundleft, pos);
            if (indexleft < 0) {
                indexleft = -indexleft - 2 +1;
            }
            else if (indexleft > -1) {
                indexleft = indexleft +1;
            }


            int indexright = Arrays.binarySearch(boundright, pos);
            if (indexright < 0){
                indexright = -indexright-1;
            }
            else if (indexright > -1){
                indexright = indexright +1;
            }


            for (int j = indexright; j < indexleft ; j++) {
                count[j]++; //每个Bin 里面的变异个数
                value[j].add(va); //每个bin里面的值的集合
            }
        }

        //计算每个Bin里面的0 1 2 NA 的个数
        String[] mean = new String[bound.length];
        for (int i = 0; i < value.length; i++) {
            mean[i] = new AoMath().getResidualHeterozygosity(value[i]);
        }

        try {
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")){
                bw = IOUtils.getTextWriter(outfileS);
            }
            if (outfileS.endsWith(".txt.gz")){
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tHETEROZYGOSITY");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(chr).append("\t").append(bound[i][0]).append("\t").append(bound[i][1]).append("\t").append(count[i])
                        .append("\t").append(mean[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println( "Bin calculation is completed at "+ outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据 chr pos 和 value 来确定,返回 每个window内的变异个数，以及pos对应值的集合的平均值
     *
     * @param chr
     * @param hm
     */
    public void calwindowstep(String chr, HashMap<Integer,String> hm, int window, int step, String outfileS){
        // 1.将pos转为为list,找到最大值，根据最大值确定bin的数目
        // 2.建立count数目和list数组，对pos进行循环，找到每个bin的左边的数目和value的集合，求这个集合的平均值，最大值，方差等
        // 3.输出，每个bin的值
        //************************
        List<Integer> posl= new ArrayList<Integer>(hm.keySet());
        Collections.sort(posl);
        int posmax = Collections.max(posl);
        int[][] bound = this.initializeWindowStep(posmax, window,step);

        int count[] = new int[bound.length]; //查看每个bin里面的变异个数
        TDoubleArrayList[] value = new TDoubleArrayList[bound.length]; //每个Bin里面的值的集合 List
        int[] boundright = new int[bound.length]; //只看左边的bound
        int[] boundleft = new int[bound.length]; //右边的bound
        for (int i = 0; i < bound.length; i++) { //每个bound的左边
            boundleft[i] = bound[i][0];
            boundright[i] = bound[i][1];
            value[i] = new TDoubleArrayList(); //对每一个bin中的List进行初始化
        }


        for (int i = 0; i < posl.size(); i++) {
            double v = Double.parseDouble(hm.get(posl.get(i)));

            int indexleft = Arrays.binarySearch(boundleft, posl.get(i));
            if (indexleft < 0) {
                indexleft = -indexleft - 2 +1;
            }
            else if (indexleft > -1) {
                indexleft = indexleft +1;
            }


            int indexright = Arrays.binarySearch(boundright, posl.get(i));
            if (indexright < 0){
                indexright = -indexright-1;
            }
            else if (indexright > -1){
                indexright = indexright +1;
            }


            for (int j = indexright; j < indexleft ; j++) {
                count[j]++; //每个Bin 里面的变异个数
                value[j].add(v); //每个bin里面的值的集合
            }
        }

        //计算每个Bin里面的平均值
        String[] mean = new String[bound.length];
        for (int i = 0; i < value.length; i++) {
            mean[i] = new AoMath().getRelativeMean(value[i]);
        }

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tHETEROZYGOSITY");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(chr).append("\t").append(bound[i][0]).append("\t").append(bound[i][1]).append("\t").append(count[i])
                        .append("\t").append(mean[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println( "It is completed at "+ outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * return the bound needed
     *
     * @param chrLength
     * @param windowSize
     * @param windowStep
     * @return
     */
    private int[][] initializeWindowStep (int chrLength, int windowSize, int windowStep) {

        TIntArrayList startList = new TIntArrayList();
        TIntArrayList endList = new TIntArrayList();
        int start = 1;
        int end = start+windowSize;
        while (start < chrLength) {
            startList.add(start);
            endList.add(end);
            start+=windowStep;
            end = start+windowSize;
        }

        int[][] bound = new int[startList.size()][2];
        for (int i = 0; i < startList.size(); i++) {
            bound[i][0] = startList.get(i);
            bound[i][1] = endList.get(i);
        }
        return bound;
    }


    /**
     * 根据 chr pos 和 value 来确定
     *
     * @param chr
     * @param hm
     */
    public void calwindow(String chr, HashMap<Integer,String> hm, int window, String outfileS){
        // 1.将pos转为为list,找到最大值，根据最大值确定bin的数目
        // 2.建立count数目和list数组，对pos进行循环，找到每个bin的左边的数目和value的集合，求这个集合的平均值，最大值，方差等
        // 3.输出，每个bin的值
        //************************

        List<Integer> posl= new ArrayList<Integer>(hm.keySet());
        Collections.sort(posl);
        int posmax = Collections.max(posl);
        int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(posmax, window);


        int count[] = new int[bound.length]; //查看每个bin里面的变异个数
        TDoubleArrayList[] value = new TDoubleArrayList[bound.length]; //每个Bin里面的值的集合 List
        int[] bounds = new int[bound.length]; //只看左边的bound
        for (int i = 0; i < bound.length; i++) { //每个bound的左边
            bounds[i] = bound[i][0];
            value[i] = new TDoubleArrayList(); //对每一个bin中的List进行初始化
        }


        for (int i = 0; i < posl.size(); i++) {
            double v = Double.parseDouble(hm.get(posl.get(i)));
            int index = Arrays.binarySearch(bounds, posl.get(i));
            if (index < 0) {
                index = -index - 2;
            }
            count[index]++; //每个Bin 里面的变异个数
            value[index].add(v); //每个bin里面的值的集合
        }

        //计算每个Bin里面的平均值
        String[] mean = new String[bound.length];
        for (int i = 0; i < value.length; i++) {
            mean[i] = new AoMath().getRelativeMean(value[i]);
        }


        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tHETEROZYGOSITY");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(chr).append("\t").append(bound[i][0]).append("\t").append(bound[i][1]).append("\t").append(count[i])
                .append("\t").append(mean[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println( "It is completed at "+ outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据数据库动态创建分组，将该分组内的所有数字建立list，进行bin的统计，并返回每个bin的比例
     *
     */
    public void getDAFtable() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014/002_basedGerpPhyloP";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014/003_basedSIFT_ratio";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014/005_basedonlyGERP";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge";

//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/004_DAFtable"; //总共的ABD
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/005_DATtable_barley_urartu";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/004_DAFtable_barley_secale_parsimony";

//        String outfileDirS = "";
        new File(outfileDirS).mkdirs();

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            try {
                //************************************ 第一阶段，定义输出输出文件，读写文件 ************************//
                String infileS = f.getAbsolutePath();
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_onlyABD.txt").getAbsolutePath(); //只能画总体的ABD六倍体
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_onlyAB.txt").getAbsolutePath(); //只能画总体的AB四倍体
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_onlyD.txt").getAbsolutePath(); //只能画总体的D二倍体

//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_Asubgenome.txt").getAbsolutePath(); //只有A亚基因组的结果
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_Bsubgenome.txt").getAbsolutePath(); //只有A亚基因组的结果
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_Dsubgenome.txt").getAbsolutePath(); //只有A亚基因组的结果

                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                //************************************ 第二阶段，创建相关变量，并读入文件判断 ************************//
                String[] group = {"Deleterious SNPs", "Nonsynonymous-tolerant SNPs", "Synonymous SNPs"};
                TDoubleArrayList[] dafABD = new TDoubleArrayList[group.length];
                for (int i = 0; i < dafABD.length; i++) { //初始化List
                    dafABD[i] = new TDoubleArrayList();
                }

                TDoubleArrayList[] dafAB = new TDoubleArrayList[group.length];
                for (int i = 0; i < dafAB.length; i++) { //初始化List
                    dafAB[i] = new TDoubleArrayList();
                }

                TDoubleArrayList[] daf = new TDoubleArrayList[group.length];
                for (int i = 0; i < daf.length; i++) { //初始化List
                    daf[i] = new TDoubleArrayList();
                }

                double sift = Double.NaN;
                double gerp = Double.NaN;
                double phylop = Double.NaN;
                String temp = null;
                String header = br.readLine();
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) { //我想得到 DAF_ABD的pos集合， DAF_AB的pos集合，以及合并数据的DAF集合。每个集合又分为3类，一类是同义突变，一类是非同义突变，sift值小于0.05，一类是非同义突变，sfft大于0.05
                    l = PStringUtils.fastSplit(temp);
                    int chrID = Integer.parseInt(l.get(1));
                    int posID = Integer.parseInt(l.get(2));
                    String chr = RefV1Utils.getChromosome(chrID,posID);
//                    if(chr.contains("D"))continue; //只能用于AABB总体画图
//                    if(chr.contains("A") || chr.contains("B"))continue; //只能用于DD总体画图

//                    if(chr.contains("B") || chr.contains("D"))continue; //只能用于A亚基因组
//                    if(chr.contains("A") || chr.contains("D"))continue; //只能用于B亚基因组
                    if(chr.contains("A") || chr.contains("B"))continue; //只能用于D亚基因组
                    System.out.println(temp);
                    String type = l.get(12);
                    String siftscore = l.get(13);
                    String gerpscore = l.get(20);
//                    String phylopscore = l.get(18);

//                    String DAF_ABD = l.get(26); //大麦黑麦为ancestral allele计算的DAF值
//                    String DAF_AB = l.get(27); //大麦黑麦为ancestral allele计算的DAF值
//                    String DAF = l.get(25); //大麦黑麦为ancestral allele计算的DAF值

//                    String DAF_ABD = l.get(29); //大麦乌拉尔图为ancestral allele计算的DAF值
//                    String DAF_AB = l.get(30); //大麦乌拉尔图为ancestral allele计算的DAF值
//                    String DAF = l.get(28); //大麦乌拉尔图为ancestral allele计算的DAF值

                    String DAF_ABD = l.get(33); //大麦黑麦Pasimony法为ancestral allele计算的DAF值
                    String DAF_AB = l.get(34); //大麦黑麦Pasimony法为ancestral allele计算的DAF值
                    String DAF = l.get(32); //大麦黑麦Pasimony法为ancestral allele计算的DAF值

                    //如果变异类型是同义突变，那么就不用做任何判断；直接加上分组 Synonymous 并写入
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值小于0.5，gerp和phylop存在，且gerp大于1，且phylop大于0.5；那么加上分组 Deleterious 并写入； gerp 值和 phylop值不满足条件的，那么就不进行分组
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值大于0.5，那么加上分组 Nonsynonymous_tolerent 并写入
                    //如果变异类型是非同义突变，SIFT值不存在，那么不分组 不写入
                    if (type.equals("SYNONYMOUS")) {
                        if (!DAF.startsWith("N")) { //说明是有值的
                            double value = Double.parseDouble(DAF);
                            daf[2].add(value);
                        }
                        if (!DAF_ABD.startsWith("N")) { //说明是有值的
                            double value = Double.parseDouble(DAF_ABD);
                            dafABD[2].add(value);
                        }
                        if (!DAF_AB.startsWith("N")) { //说明是有值的
                            double value = Double.parseDouble(DAF_AB);
                            dafAB[2].add(value);
                        }

                    }
                    if (type.equals("NONSYNONYMOUS")) {
                        if (!siftscore.startsWith("N")) {
                            sift = Double.parseDouble(siftscore);
                            if (sift < 0.05) {
                                //添加gerp phyloP分组信息
//                                if (!gerpscore.startsWith("N") && (!phylopscore.startsWith("N"))) { //均有值存在
//                                    gerp = Double.parseDouble(gerpscore);
//                                    phylop = Double.parseDouble(phylopscore);
//                                    if (gerp > 1 && (phylop > 0.5)) {
//                                        if (!DAF.startsWith("N")) { //说明是有值的
//                                            double value = Double.parseDouble(DAF);
//                                            daf[0].add(value);
//                                        }
//                                        if (!DAF_ABD.startsWith("N")) { //说明是有值的
//                                            double value = Double.parseDouble(DAF_ABD);
//                                            dafABD[0].add(value);
//                                        }
//                                        if (!DAF_AB.startsWith("N")) { //说明是有值的
//                                            double value = Double.parseDouble(DAF_AB);
//                                            dafAB[0].add(value);
//                                        }
//                                    }
//                                }


                                //不添加gerp phyloP分组信息
//**************************** 可供选择 *********************************************** //
//                                if (!DAF.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF);
//                                    daf[0].add(value);
//                                }
//                                if (!DAF_ABD.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_ABD);
//                                    dafABD[0].add(value);
//                                }
//                                if (!DAF_AB.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_AB);
//                                    dafAB[0].add(value);

//**************************** 可供选择 *********************************************** //

                                if (!gerpscore.startsWith("N")) { //均有值存在
                                    gerp = Double.parseDouble(gerpscore);
                                    if (gerp > 1) {
                                        if (!DAF.startsWith("N")) { //说明是有值的
                                            double value = Double.parseDouble(DAF);
                                            daf[0].add(value);
                                        }
                                        if (!DAF_ABD.startsWith("N")) { //说明是有值的
                                            double value = Double.parseDouble(DAF_ABD);
                                            dafABD[0].add(value);
                                        }
                                        if (!DAF_AB.startsWith("N")) { //说明是有值的
                                            double value = Double.parseDouble(DAF_AB);
                                            dafAB[0].add(value);
                                        }
                                    }
                                }

                            } else { //sift值大于等于0.05

                                if (!DAF.startsWith("N")) { //说明是有值的
                                    double value = Double.parseDouble(DAF);
                                    daf[1].add(value);
                                }
                                if (!DAF_ABD.startsWith("N")) { //说明是有值的
                                    double value = Double.parseDouble(DAF_ABD);
                                    dafABD[1].add(value);
                                }
                                if (!DAF_AB.startsWith("N")) { //说明是有值的
                                    double value = Double.parseDouble(DAF_AB);
                                    dafAB[1].add(value);
                                }
                            }
                        }
                    }
                }
                br.close(); //文件阅读完毕！！

                //************************************ 第三阶段，开始返回bin的比例 ************************//

                List[] l1 = this.mkBarplotofDAF(dafABD[0], 20);
                List[] l2 = this.mkBarplotofDAF(dafABD[1], 20);
                List[] l3 = this.mkBarplotofDAF(dafABD[2], 20);
                List[] l4 = this.mkBarplotofDAF(dafAB[0], 20);
                List[] l5 = this.mkBarplotofDAF(dafAB[1], 20);
                List[] l6 = this.mkBarplotofDAF(dafAB[2], 20);
                List[] l7 = this.mkBarplotofDAF(daf[0], 20);
                List[] l8 = this.mkBarplotofDAF(daf[1], 20);
                List[] l9 = this.mkBarplotofDAF(daf[2], 20);
                //************************************ 第四阶段，开始写出文件 ************************//
                bw.write("Xaxes\tDAF_ABD\tDAF_AB\tDAF\tGroup");
                bw.newLine();
                for (int i = 0; i < l1[0].size(); i++) {
                    bw.write(String.valueOf(l1[0].get(i)) + "\t" + String.valueOf(l1[1].get(i)) + "\t" + String.valueOf(l4[1].get(i))+ "\t" + String.valueOf(l7[1].get(i)) + "\t" + group[0]);
                    bw.newLine();
                }
                for (int i = 0; i < l1[0].size(); i++) {
                    bw.write(String.valueOf(l1[0].get(i))+ "\t" + String.valueOf(l2[1].get(i)) + "\t" + String.valueOf(l5[1].get(i))+ "\t" + String.valueOf(l8[1].get(i))+ "\t" + group[1]);
                    bw.newLine();
                }
                for (int i = 0; i < l1[0].size(); i++) {
                    bw.write(String.valueOf(l1[0].get(i)) + "\t"+ String.valueOf(l3[1].get(i)) + "\t" + String.valueOf(l6[1].get(i))+ "\t" + String.valueOf(l9[1].get(i))+ "\t" + group[2]);
                    bw.newLine();
                }
                

//                bw.write("Xaxes\tDeleterious_SNPs\tNonsynonymous_tolerant_SNPs\tSynonymous_SNPs\tGroup");
//                bw.newLine();
//                for (int i = 0; i < l1[0].size(); i++) {
//                    bw.write(String.valueOf(l1[0].get(i)) + "\t" + String.valueOf(l1[1].get(i)) + "\t" + String.valueOf(l2[1].get(i))+ "\t" + String.valueOf(l3[1].get(i)) + "\t" + "Hexaploid");
//                    bw.newLine();
//                }
//                for (int i = 0; i < l1[0].size(); i++) {
//                    bw.write(String.valueOf(l1[0].get(i))+ "\t" + String.valueOf(l4[1].get(i)) + "\t" + String.valueOf(l5[1].get(i))+ "\t" + String.valueOf(l6[1].get(i))+ "\t" + "Tetraploid");
//                    bw.newLine();
//                }
//                for (int i = 0; i < l1[0].size(); i++) {
//                    bw.write(String.valueOf(l1[0].get(i)) + "\t"+ String.valueOf(l7[1].get(i)) + "\t" + String.valueOf(l8[1].get(i))+ "\t" + String.valueOf(l9[1].get(i))+ "\t" + "Diploid");
//                    bw.newLine();
//                }
                
                bw.flush();
                bw.close();

                //
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    public List[] mkBarplotofDAF(TDoubleArrayList dafList, int bins) {
        HashMap<Double, Double> hm = new HashMap<>();
        List[] l = new List[2];
        for (int i = 0; i < l.length; i++) {
            l[i] = new ArrayList();
        }
//        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble("1");

        //先建立bound数组,假如 bin=20,那么每个bin的大小是0.05 bound[0] =0; bound[1]=0.05, bound[19]=0.95 即只取了区间的左边开始位置
        double[] bound = new double[bins];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double) length / bins * i;
        }

        double[] daf = new double[bins]; //确定落入每个区间的个数，进行计算，20个bin，20个count数字
        for (int i = 0; i < dafList.size(); i++) { //对每个daf值进行区间的判断
            double value = dafList.get(i);
            int index = Arrays.binarySearch(bound, value);
            if (index < 0) {
                index = -index - 2;
            }
            //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
            //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
            daf[index]++; //值落入第i种变异的第index个区间的个数
        }
        //开始计算每个区间落入点的比例
        for (int i = 0; i < daf.length; i++) {
            daf[i] = daf[i] / dafList.size();
        }
        System.out.println("This list is " + dafList.size() + "  size");
        //开始写出文件
        try {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < bound.length; i++) {
                String coordinate = String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2);
                String proportion = String.format("%.4f", daf[i]);
                l[0].add(coordinate);
                l[1].add(proportion);
                System.out.println(coordinate + "   " + proportion);
            }
            
            
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return l;
    }

    /**
     *
     * @param infileS
     * @param outfileDirS
     * @param binNum the number of bins that would be divided
     */
    public void mkBarplotofMAF_single(String infileS, String outfileDirS, String binNum, String max) {
        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble(max);
        new File(outfileDirS).mkdirs();
        String outfileS = null;
        BufferedReader br = null;
        if (infileS.endsWith(".txt")) {
            br = IOUtils.getTextReader(infileS);
            outfileS = new File(outfileDirS, new File(infileS).getName().replaceFirst("txt", bins + "bins" + ".Table.txt")).getAbsolutePath();
        } else if (infileS.endsWith(".txt.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
            outfileS = new File(outfileDirS, new File(infileS).getName().replaceFirst("txt.gz", bins + "bins" + ".Table.txt")).getAbsolutePath();
        }
        //先建立bound数组
        double[] bound = new double[bins];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double) length / bins * i;
        }

        double[] maf = new double[bins];
        TDoubleArrayList mafList = new TDoubleArrayList();
        RowTable<String> t = new RowTable<>(infileS);
        int count = t.getRowNumber();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCell(i, 1).equals("NA")) { //MAF值所在的那一列
                continue;
            }
            double value = t.getCellAsDouble(i, 1); //MAF值所在的那一列
            mafList.add(value);
            int index = Arrays.binarySearch(bound, value);
            if (index < 0) {
                index = -index - 2;
            }
            //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
            //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
            maf[index]++; //值落入第i种变异的第index个区间的个数
        }
        //开始计算每个区间落入点的比例
        for (int i = 0; i < maf.length; i++) {
            maf[i] = maf[i] / mafList.size();
        }
        System.out.println(mafList.size() + "  size");
        //开始写出文件
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Maf\tDensity");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("%.2f", bound[i] + length / bins / 2)).append("\t").append(String.format("%.4f", maf[i]));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(new File(infileS).getName() + " is completed at " + outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     *
     * @param infileDirS
     * @param outfileDirS
     * @param binNum the number of bins that would be divided
     */
    public void mkBarplotofMAF(String infileDirS, String outfileDirS, String binNum, String max) {
        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble(max);
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt", bins + "bins" + ".Table.txt")).getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt.gz", bins + "bins" + ".Table.txt")).getAbsolutePath();
            }
            //先建立bound数组
            double[] bound = new double[bins];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double) length / bins * i;
            }

            double[] maf = new double[bins];
            TDoubleArrayList mafList = new TDoubleArrayList();
            RowTable<String> t = new RowTable<>(infileS);
            int count = t.getRowNumber();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 7).equals("NA")) { //MAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 7); //MAF值所在的那一列
                mafList.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                maf[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < maf.length; i++) {
                maf[i] = maf[i] / mafList.size();
            }
            System.out.println(mafList.size() + "  size");
            //开始写出文件
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Daf\tDensity");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2)).append("\t").append(String.format("%.4f", maf[i]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * make bin table, 一个文件制作一个Binfile，有起始和终止位置,这里总长度根据文件的pos信息获得
     *
     * @param infileDirS
     * @param outfileDirS
     * @param binsize
     */
    public void mkBintable(String infileDirS, String outfileDirS, String binsize) {
        int binSize = Integer.parseInt(binsize);
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "." + binSize / 1000000 + "M" + ".binTable.txt")).getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt.gz", "." + binSize / 1000000 + "M" + ".binTable.txt")).getAbsolutePath();
            }

            RowTable<String> t = new RowTable(infileS);
            //染色体的长度是最后一行pos的位置
            int chrlength = Integer.valueOf(t.getCell(t.getRowNumber() - 1, 1)); //t.getRowNumber()是文件的行数，不包括header。 这里getCell得到的是索引
            String Chr = t.getCell(0, 0);

            int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(chrlength, binSize);
            int count[] = new int[bound.length];
            int[] bounds = new int[bound.length];
            for (int i = 0; i < bound.length; i++) { //每个bound的左边
                bounds[i] = bound[i][0];
            }
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(bounds, Integer.valueOf(t.getCell(i, 1)));
                if (index < 0) {
                    index = -index - 2;
                }
                count[index]++;
            }

            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("CHROM\tBIN_START\tSNP_COUNT\tVARIANTS.KB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    double variant = (double) count[i] / (double) 1000;
                    //sb.append(Chr).append("\t").append(bound[i][0]).append("\t").append(bound[i][1]).append("\t").append(variant);
                    sb.append(Chr).append("\t").append(bound[i][0]).append("\t").append(count[i]).append("\t").append(variant);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName() + " is completed.");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * extract pos info from vcf file. eg: vcf ---- Chr Pos
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mkHapPos(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".pos.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".pos.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnt = 0;
                bw.write("Chr\tPos\n");
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temp = temp.substring(0, 40); //肯定够                 
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt % 1000 == 0) {
                        System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    }
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

}
