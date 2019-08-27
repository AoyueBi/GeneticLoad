/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class FilterVCF {
    
    public FilterVCF(){
        //this.statVcfCoverage();
        //this.subsetCovevsSDvsPV();
        //this.statVcfPValue();
        this.subsetCovevsSDvsPV();
        
    }
    
    
    /**
     * 对已经生成的2万行结果进行随机抽样，抽出5000行进行接下来的随机分析
     */
    public void subsetCovevsSDvsPV() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/abd/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/abd/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/ab/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/ab/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/d/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/d/";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/mergeTaxa/";


        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".depth.txt.gz", ".depth_5000sites.txt")).getAbsolutePath();
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int totallines = t.getRowNumber();
            double ratio = (double) 5000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数
            for (int i = 0; i < t.getRowNumber(); i++) {
                double r = Math.random();
                if (r > ratio) {

                } else {
                    ifOut[i] = true;
                }
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
        });
    }
    
    /**
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计和Pvalue获取，制成一个表格，chr pos averageDepth SD PValue
     * IfPVzero taxa1Depth ..... 注意表格不要以#开头，否则被注释，看不到
     * 代码，即在之前depth和SD的基础上又加了2列，第一列是PValue 第二列是PValue是否是0的判断,如果是，那么值为1，如果不是那么值为0
     */
    public void statVcfPValue() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/abd/001_subset/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/abd/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/ab/001_subsetVCF/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/ab/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/d/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";


        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);

        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            //infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/chr005.merge.vcf.gz";
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".depth.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + "AverageDepth\tSD\tPValue\tIfPVzero");
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                    bw.write("\t" + taxa[i]);
                }
                bw.newLine();
                //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();

                    cnt++; // 对snp开始计数
                    if (cnt % 1000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    String pvalue = l.get(7).split("PV=")[1].split(";")[0];
                    String ifPVzero = null;
                    if (pvalue.equals("0.0")) {
                        ifPVzero = "1";
                    } else {
                        ifPVzero = "0";
                    }
                    for (int i = 0; i < taxa.length; i++) {

                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    //计算完毕，接下来开始写入文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(String.format("%.3f", relativeMean)).append("\t").append(String.format("%.3f", sd))
                            .append("\t").append(pvalue).append("\t").append(ifPVzero);
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(String.format("%.0f", depthList.get(i)));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }
    
    /**
     * 对已经生成的2万行结果进行随机抽样，抽出5000行进行接下来的随机分析
     */
    public void subsetCoveVSSD() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/abd/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/002_depthCalSample/abd/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/ab/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/002_depthCalSample/ab/";
        
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/d/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/002_depthCalSample/d/";


        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".depth.txt.gz", ".depth_5000sites.txt")).getAbsolutePath();
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int totallines = t.getRowNumber();
            double ratio = (double) 5000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数
            for (int i = 0; i < t.getRowNumber(); i++) {
                double r = Math.random();
                if (r > ratio) {

                } else {
                    ifOut[i] = true;
                }
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
        });
    }
    
    
    /**
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计，制成一个表格，chr pos averageDepth SD taxa1Depth
     * ..... 注意表格不要以#开头，否则被注释，看不到
     */
    public void statVcfCoverage() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/abd/001_subset/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/abd/";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/ab/001_subsetVCF/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/ab/";
        
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/d/";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".depth.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + "AverageDepth\tSD");
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                    bw.write("\t" + taxa[i]);
                }
                bw.newLine();
                //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    cnt++; // 对snp开始计数
                    if (cnt % 1000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    for (int i = 0; i < taxa.length; i++) {

                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(relativeMean).append("\t").append(sd);
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(String.format("%.0f", depthList.get(i)));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }
    
}
