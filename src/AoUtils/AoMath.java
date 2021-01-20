/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class AoMath {

    public AoMath() {

    }


    /**
     * Goal: 非连续抽样，从数组sampleArray中，非连续抽sampleSize个元素，（10，30，80..n）
     * 返回一个list类型的数组，数组元素1包含1个抽样值，元素2包含2个抽样值，元素三包含3个抽样值，以此类推。
     * @param sampleSizeArray
     * @param sampleArray
     * @return
     */
    public static List<String>[] noncontinuousRandom(String[] sampleArray, TIntArrayList sampleSizeArray){

        List<String>[] out = new List[sampleSizeArray.size()];
        for (int i = 0; i < out.length ; i++) {
            out[i] = new ArrayList<>();
        }

        for (int i = 0; i < sampleSizeArray.size(); i++) { //抽10个 20个
            int currentsize = sampleSizeArray.get(i);
            Set<String> set = new HashSet<>(currentsize);

            for (int j = 0; j < 100000; j++) {
                Random r=new Random();
                String element = sampleArray[r.nextInt(sampleArray.length)];
                set.add(element);
                if (set.size()>=currentsize)break;
            }

            List<String> l = new ArrayList<>(set);
            Collections.sort(l);
            out[i] = l;

            for (int j = 0; j < l.size(); j++) {
                System.out.print(l.get(j) + "\t");
            }
            System.out.println();
        }

        return out;
    }

    /**
     * Goal: 非连续抽样，从数组sampleArray中，非连续抽sampleSize个元素，（10，30，80..n）
     * 返回一个list类型的数组，数组元素1包含1个抽样值，元素2包含2个抽样值，元素三包含3个抽样值，以此类推。
     * @param sampleSizeArray
     * @param sampleArray
     * @return
     */
    public static List<String>[] noncontinuousRandom(String[] sampleArray,int[] sampleSizeArray){

        List<String>[] out = new List[sampleSizeArray.length];
        for (int i = 0; i < out.length ; i++) {
            out[i] = new ArrayList<>();
        }

        for (int i = 0; i < sampleSizeArray.length; i++) { //抽10个 20个
            int currentsize = sampleSizeArray[i];
            Set<String> set = new HashSet<>(currentsize);

            for (int j = 0; j < 100000; j++) {
                Random r=new Random();
                String element = sampleArray[r.nextInt(sampleArray.length)];
                set.add(element);
                if (set.size()>=currentsize)break;
            }

            List<String> l = new ArrayList<>(set);
            Collections.sort(l);
            out[i] = l;

            for (int j = 0; j < l.size(); j++) {
                System.out.print(l.get(j) + "\t");
            }
            System.out.println();
        }

        return out;
    }

    /**
     * Goal: 连续抽样，从数组sampleArray中，连续抽sampleSize个元素，（1..n）
     * 返回一个list类型的数组，数组元素1包含1个抽样值，元素2包含2个抽样值，元素三包含3个抽样值，以此类推。
     * @param sampleSize
     * @param sampleArray
     * @return
     */
    public static List<String>[] continuousRandom(String[] sampleArray,int sampleSize){

        List<String>[] out = new List[sampleSize];
        for (int i = 0; i < out.length ; i++) {
            out[i] = new ArrayList<>();
        }

        for (int i = 0; i < sampleSize; i++) { //抽1个 2个 3个 4个 5个 6个。。。
            int currentsize = i+1;
            Set<String> set = new HashSet<>(currentsize);

            for (int j = 0; j < 100000; j++) {
                Random r=new Random();
                String element = sampleArray[r.nextInt(sampleArray.length)];
                set.add(element);
                if (set.size()>=currentsize)break;
            }

            List<String> l = new ArrayList<>(set);
            Collections.sort(l);
            out[i] = l;

            for (int j = 0; j < l.size(); j++) {
                System.out.print(l.get(j) + "\t");
            }
            System.out.println();
        }

        return out;
    }

    /**
     * Goal: 非连续抽样，从数组sampleArray中，非连续抽sampleSize个元素，（10，30，80..n）
     * 返回一个list类型的数组，数组元素1包含1个抽样值，元素2包含2个抽样值，元素三包含3个抽样值，以此类推。
     * @param sampleSizeArray
     * @param sampleArray
     * @return
     */
    public static List<String>[] noncontinuousRandom_withoutReplacement(String[] sampleArray, TIntArrayList sampleSizeArray){

        List<String>[] out = new List[sampleSizeArray.size()];
        for (int i = 0; i < out.length ; i++) {
            out[i] = new ArrayList<>();
        }

        for (int i = 0; i < sampleSizeArray.size(); i++) { //抽10个 20个
            int currentsize = sampleSizeArray.get(i);
            Set<String> set = new HashSet<>(currentsize);

            ////如果i>1。那么抽样的时候就将上一抽样的taxa全部加入下一抽样，剩下的随机抽样补齐。这样就解决了不放回的问题
            if (currentsize>1){
                for (int j = 0; j < out[i-1].size(); j++) {
                    set.add(out[i-1].get(j));
                }
            }

            for (int j = 0; j < 100000; j++) {
                Random r=new Random();
                String element = sampleArray[r.nextInt(sampleArray.length)];
                set.add(element);
                if (set.size()>=currentsize)break;
            }

            List<String> l = new ArrayList<>(set);
            Collections.sort(l);
            out[i] = l;

            for (int j = 0; j < l.size(); j++) {
                System.out.print(l.get(j) + "\t");
            }
            System.out.println();
        }

        return out;
    }


    /**
     * Goal: 非连续抽样，从数组sampleArray中，非连续抽sampleSize个元素，（10，30，80..n）
     * 返回一个list类型的数组，数组元素1包含1个抽样值，元素2包含2个抽样值，元素三包含3个抽样值，以此类推。
     * @param sampleSizeArray
     * @param sampleArray
     * @return
     */
    public static List<String>[] noncontinuousRandom_withoutReplacement(String[] sampleArray,int[] sampleSizeArray){

        List<String>[] out = new List[sampleSizeArray.length];
        for (int i = 0; i < out.length ; i++) {
            out[i] = new ArrayList<>();
        }

        for (int i = 0; i < sampleSizeArray.length; i++) { //抽10个 20个
            int currentsize = sampleSizeArray[i];
            Set<String> set = new HashSet<>(currentsize);

            ////如果i>1。那么抽样的时候就将上一抽样的taxa全部加入下一抽样，剩下的随机抽样补齐。这样就解决了不放回的问题
            if (currentsize>1){
                for (int j = 0; j < out[i-1].size(); j++) {
                    set.add(out[i-1].get(j));
                }
            }

            for (int j = 0; j < 100000; j++) {
                Random r=new Random();
                String element = sampleArray[r.nextInt(sampleArray.length)];
                set.add(element);
                if (set.size()>=currentsize)break;
            }

            List<String> l = new ArrayList<>(set);
            Collections.sort(l);
            out[i] = l;

            for (int j = 0; j < l.size(); j++) {
                System.out.print(l.get(j) + "\t");
            }
            System.out.println();
        }

        return out;
    }


    /**
     * Goal: 连续不放回抽样，从数组sampleArray中，连续抽sampleSize个元素，（1...n）
     * Note: 连续抽取10次，注意
     * 返回一个list类型的数组，数组元素1包含1个抽样值，元素2包含2个抽样值，元素三包含3个抽样值，以此类推。
     * @param sampleSize
     * @param sampleArray
     * @return
     */
    public static List<String>[] continuousRandom_withoutReplacement(String[] sampleArray,int sampleSize){

        List<String>[] out = new List[sampleSize];
        for (int i = 0; i < out.length ; i++) {
            out[i] = new ArrayList<>();
        }

        for (int i = 0; i < sampleSize; i++) { //抽1个 2个 3个 4个 5个 6个。。。
            int currentsize = i+1;
            Set<String> set = new HashSet<>(currentsize);
            ////如果i>1。那么抽样的时候就将上一抽样的taxa全部加入下一抽样，剩下的随机抽样补齐。这样就解决了不放回的问题
            if (currentsize>1){
                for (int j = 0; j < out[i-1].size(); j++) {
                    set.add(out[i-1].get(j));
                }
            }

            for (int j = 0; j < 100000; j++) {
                Random r=new Random();
                String element = sampleArray[r.nextInt(sampleArray.length)];
                set.add(element);
                if (set.size()>=currentsize)break;
            }

            List<String> l = new ArrayList<>(set);
            Collections.sort(l);
            out[i] = l;

            for (int j = 0; j < l.size(); j++) {
                System.out.print(l.get(j) + "\t");
            }
            System.out.println();
        }

        return out;
    }



    public static double[] NormalizeScore(double[] values){
        double[] out = new double[values.length];
        TDoubleArrayList valuez = new TDoubleArrayList();

        double max = StatUtils.max(values);
        for ( double value: values ) {
            double stdscore = (double) value*50/max;
            valuez.add(stdscore);
        }
        out = valuez.toArray();

        return out;
    }

    public static double[] ZScore(double[] values){
        double[] out = new double[values.length];
        TDoubleArrayList valuez = new TDoubleArrayList();
        double variance = StatUtils.populationVariance(values);
        double sd = Math.sqrt(variance);
        double mean = StatUtils.mean(values);
//        NormalDistribution nd = new NormalDistribution();
        for ( double value: values ) {
            double stdscore = (value-mean)/sd;
            valuez.add(stdscore);
//            double sf = 1.0 - nd.cumulativeProbability(Math.abs(stdscore)); //计算概率值
//            System.out.println("" + stdscore + " " + sf);
        }
        out = valuez.toArray();

        return out;
    }

    /**
     * 计算一组数据的z-score和p value
     */
    private void run() {
        double[] values = {9967,11281,10752,10576,2366,11882,11798};
        double variance = StatUtils.populationVariance(values);
        double sd = Math.sqrt(variance);
        double mean = StatUtils.mean(values);
        NormalDistribution nd = new NormalDistribution();
        for ( double value: values ) {
            double stdscore = (value-mean)/sd;
            double sf = 1.0 - nd.cumulativeProbability(Math.abs(stdscore));
            System.out.println("" + stdscore + " " + sf);
        }
    }

    /**
     * 获取一个表格中，第1列分类下的，第2列的每个分类的set
     * 例如：一个表格有10行
     * 第一列包含3个因子：AABBDD AABB DD
     * 第4列
     */
    public void getnlevelsforEachGroup(String infileS, int group1ColumnIndex, int group2ColumnIndex){
        //第一步：获取第一列的set集合
        List<String> lset1 = AoFile.getStringListbySet(infileS,group1ColumnIndex); //代表set转化而成的集合
        List<String> lset2 = AoFile.getStringListbySet(infileS,group2ColumnIndex); //代表set转化而成的集合
        System.out.println("******** The " + group1ColumnIndex + " column with factor " + lset1.size());
        System.out.println("******** The " + group2ColumnIndex + " column with factor " + lset2.size());
        //第二步：打印每一列含有的因子
        for (int i = 0; i < lset1.size(); i++) {
            System.out.print(lset1.get(i) + "\t");
        }
        System.out.println("");
        for (int i = 0; i < lset2.size(); i++) {
            System.out.print(lset2.get(i) + "\t");
        }
        System.out.println("");
        //第三步：获取第一列每个因子包含的集合

        int[][] sum = new int[lset2.size()][lset1.size()];
        try{
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String element1 = l.get(group1ColumnIndex);
                String element2 = l.get(group2ColumnIndex);
                int index1 = Collections.binarySearch(lset1,element1);
                int index2 = Collections.binarySearch(lset2,element2);
                sum[index2][index1]++;
            }
            br.close();

            //写出文件
            System.out.println("*******************************************************");
            System.out.println("************** start to write the table ***************");
            System.out.println("*******************************************************");

            //以下是横表，比较长，不实用。
//            System.out.print("sum");
//            for (int i = 0; i < lset2.size(); i++) {
//                System.out.print("\t" + lset2.get(i));
//            }
//            System.out.println("\tNlevel");
//
//            //统计每行的不是0的因子个数
//            for (int i = 0; i < lset1.size(); i++) {
//                System.out.print(lset1.get(i));
//                int numfactor = 0;
//                for (int j = 0; j < lset2.size(); j++) {
//                    System.out.print("\t" + sum[i][j]);
//                    if (sum[i][j] !=0) numfactor++;
//                }
//                System.out.println("\t" + numfactor);
//            }

            //以下是长表，比较短，实用。
            int[] numfactor = new int[lset1.size()];
            System.out.print("region");
            for (int i = 0; i < lset1.size(); i++) {
                System.out.print("\t" + lset1.get(i));
            }
            System.out.println("");

            for (int i = 0; i < lset2.size(); i++) {
                System.out.print(lset2.get(i));
                for (int j = 0; j < lset1.size(); j++) {
                    System.out.print("\t" + sum[i][j]);
                    if (sum[i][j] !=0) numfactor[j]++;
                }
                System.out.println("");
            }
            System.out.print("Nlevel\t");
            for (int i = 0; i < numfactor.length; i++) {
                System.out.print(numfactor[i] + "\t");
            }
            System.out.println("");
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 将注释库中DAF是0或者1的位点一环为NA,NA的位点依旧是NA
     */
    public static String replace10toNA(String value){
        String out = null;
        double daf = Double.NaN;
        if (value.startsWith("N")){
            out = "NA";
        }
        if(!value.startsWith("N")){
            daf = Double.parseDouble(value);
            if (daf == 1 || daf== 0){
                out = "NA";
            }
            else{
                out = value;
            }
        }
        return out;
    }

    /**
     * 求所有集合的总和
     * @param l
     * @return
     */
    public static int listSum_byint(List<Integer> l){
        int out = 0;
        for (int i = 0; i < l.size() ; i++) {
            out = out + l.get(i);
        }

        return out;
    }

    /**
     *
     * @return
     */
    public static List<Double> topK(){
        //我想根据pos xpclr 得到
        List<Double> out = new ArrayList<>();

        List<Double> l = new ArrayList<Double>();
//        double[] array = {1,5,2,9,4,7,5.6,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY};
        Double[] array = {1.0, 5.0, 2.0, 9.0, 4.0, 7.0, 5.6, 8.8, 9.0, 3.2, 4.5, 6.6, 1.2, 9.4, 4.5, 6.6, 1.2, 9.4, 9.8, 7.3 };
        double k = 0.2;
        int length = array.length;
        double kk = k*length;
        l = Arrays.asList(array);
        Collections.sort(l,Collections.reverseOrder());
        for (int i = 0; i < kk ; i++) {
//            System.out.println(l.get(i));
            System.out.println(l.indexOf(l.get(i)));
        }

        return out;
    }


    /**
     * there must be no NA value
     * @param value
     * @return
     */
    public static String getRelativeMean(TDoubleArrayList value){
        String out = null;
        TDoubleArrayList list = new TDoubleArrayList();


        for (int i = 0; i < value.size(); i++) {
            if (!Double.isNaN(value.get(i))){
                list.add(value.get(i));
            }
        }

        double[] array = list.toArray();
        DescriptiveStatistics d = new DescriptiveStatistics(array);
        double relativeMean = d.getMean(); //平均值
        out = String.format("%.4f", relativeMean);
        return out;
    }


    public static String getStandardDeviation(TDoubleArrayList value){
        String out = null;
        double[] array = value.toArray();
        DescriptiveStatistics d = new DescriptiveStatistics(array);
        double sd = d.getStandardDeviation(); //标准偏差
        out = String.format("%.4f", sd);
        return out;
    }


    public String descriptiveStatistics(TDoubleArrayList value){
        String out = null;
        double[] array = value.toArray();
        DescriptiveStatistics d = new DescriptiveStatistics(array);
        double relativeMean = d.getMean(); //平均值
        double sd = d.getStandardDeviation(); //标准偏差
        double median = d.getPercentile(50); //中位数
        double x = d.getKurtosis();
        out = String.format("%.4f", relativeMean);
        return out;
    }


    /**
     * return the RH value from 0(0/0) 1(0/1) 2(1/1) NA(./.)
     *
     * @param l
     * @return
     */
    public String getResidualHeterozygosity(List<String> l){
        String out = null;
        int cnt0 = 0;
        int cnt1 = 0;
        int cnt2 = 0;
        int cntNA = 0;

        for (int i = 0; i < l.size(); i++) {
            String v = l.get(i);
            if (v.equals("0")){
                cnt0++;
            }
            else if (v.equals("1")){
                cnt1++;
            }
            else if (v.equals("2")){
                cnt2++;
            }
            else if (v.equals("NA")){
                cntNA++;
            }
        }
        int cnt = cnt0 + cnt1 + cnt2;
        Double a = (double) cnt1/cnt;
        out = String.format("%.4f",a);

        return out;
    }

    /**
     *  将group的名字用 index 数字代替
     * @param infileS
     * @param columnIndex
     * @return
     */
    public HashMap<String,Integer> setGrouptoNumber (String infileS, int columnIndex){
        HashMap<String,Integer> hm = new HashMap<>();
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(columnIndex);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "\t" + Collections.frequency(l, a));
        }
        List<String> group  = new ArrayList<>(s);
        for (int i = 0; i < group.size();i++) {
            String name = group.get(i);
            hm.put(name,i+1);
        }

        return hm;
    }

    public static void countCaseInGroup(String infileS, int columnIndex){
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(columnIndex);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "\t" + Collections.frequency(l, a));
        }
    }

    public static void countCase_fromList(List<String> l){
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "\t" + Collections.frequency(l, a));
        }
    }

    public static File countCase_fromList_outFile(List<String> l){
        File out = new File("/Users/Aoyue/Documents/countCase.txt");
        BufferedWriter bw = AoFile.writeFile(out.getAbsolutePath());
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);

        try{
            bw.write("CHROM\tCount\tSub");
            bw.newLine();
            for(String a : s){
                System.out.println(a + "\t" + Collections.frequency(l, a));
                String sub = a.substring(1);
                bw.write(a+"\t"+Collections.frequency(l, a)+"\t"+sub);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }


    /**
     * 过滤DAF_ABD等于0或者1的位点
     *
     */
    public void filterValue(String infileDirS, String outfileDirS) {
//        new AoFile().readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/010_genicSNPAnnotation_addGERPandPhyloP/chr001_SNP_anno.txt.gz");

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/004_merge";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/005_filterDAF_ABD";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach((File f) -> {

            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_filterGERPandPhyloP.txt")).getAbsolutePath();

                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                BufferedWriter bw = null; // IOUtils.getTextGzipWriter(outfileS);
                if (outfileS.endsWith(".txt")) {
                    bw = IOUtils.getTextWriter(outfileS);
                } else if (outfileS.endsWith(".txt.gz")) {
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }
                String temp = null;
                String header = br.readLine(); //读表头
                bw.write(header);
                bw.newLine();
                int cnt = 0;
                List<String> l = new ArrayList();
                String goalValue1 = null;
                String goalValue2 = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    cnt++;
//ID	Chr	Pos	Ref	Alt	Major	Minor	Maf	AAF_ABD	AAF_AB	Transcript	Region	Variant_type	SIFT_score	Ancestral	DAF	DAF_ABD	DAF_AB	Gerp	PhyloP
                    goalValue1 = l.get(18); //此处需要修改，目标值
                    goalValue2 = l.get(19);
                    if (goalValue1.startsWith("N")) {
                        continue;
                    }
                    if (goalValue2.startsWith("N")) {
                        continue;
                    }

                    double value1 = Double.parseDouble(goalValue1);
                    double value2 = Double.parseDouble(goalValue2);
//                    if (value1 == 0 || (value1 == 1)) {
//                        continue;
//                    }
//                    if (value1 < 1 || (value2 < 0.5)) {
//                        continue;
//                    }
                    if (value1 < 0.05) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();

                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
        /*==================================== 测试用 =============================================*/
    }

}
