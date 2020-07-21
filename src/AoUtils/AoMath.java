/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import gnu.trove.list.array.TDoubleArrayList;
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


    public static String getRelativeMean(TDoubleArrayList value){
        String out = null;
        double[] array = value.toArray();
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
