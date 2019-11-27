/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.Bin;
import AoUtils.CountSites;
import format.genomeAnnotation.GeneFeature;
import format.range.Range;
import format.table.RowTable;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class VariantsSum {

    public VariantsSum() {
//        this.mkSNPsummary("/data4/home/aoyue/vmap2/genotype/mergedVCF/002_biMAF0.005VCF/", "/data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/");
        //this.addAncestralAllele("/Users/Aoyue/Documents/out", "/Users/Aoyue/Documents/out1", "/Users/Aoyue/Documents/out2");
        //this.scriptAddAncAllele();
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/delSNP");
//new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/nonsyTolerantSNP");
//new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/synSNP");
        // new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/004_merge/delSNP");
        // new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/004_merge/nonsyTolerantSNP");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/004_merge/synSNP");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/delSNP/chrA.subgenome.delSNP.changeChrPos.txt");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/delSNP/chrB.subgenome.delSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/nonsyTolerantSNP/chrA.subgenome.nonsyTolerantSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/nonsyTolerantSNP/chrB.subgenome.nonsyTolerantSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/synSNP/chrA.subgenome.synSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/synSNP/chrB.subgenome.synSNP.changeChrPos.txt");
//this.mkBarplotOfSNPs();

        //this.classifySNPs("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/005_addAncestralAllele/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/");
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/delSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/delSNP");
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/nonsyTolerantSNP");
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/synSNP");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/delSNP/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/delSNP/");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/nonsyTolerantSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/nonsyTolerantSNP");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/synSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/synSNP");
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/delSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/delSNP/chrD.subgenome.delSNP.changeChrPos.txt");
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/nonsyTolerantSNP/chrD.subgenome.nonsyTolerantSNP.changeChrPos.txt");
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/synSNP/chrD.subgenome.synSNP.changeChrPos.txt");
//        this.mkBarplotOfSNPs();
//        this.mkSNPsummary_step1("/Users/Aoyue/Documents/chr001_vmap2.1_line50.vcf.gz", "/Users/Aoyue/Documents/out/chr001_vmap2.1_line50.annotation.txt.gz");
//        this.mkSNPsummary_step2("/Users/Aoyue/Documents/chr001_vmap2.1_AnnoDB_10000lines.txt", "/Users/Aoyue/Documents/chr001.wheat.ancestralAllele_1000000.txt", "/Users/Aoyue/Documents/out/chr001_vmap2.1_.txt.gz");
//    this.getCDSannotation("/Users/Aoyue/Documents/test", "/Users/Aoyue/Documents/out");
//        this.classifySNP_byPop("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/102_cdsAnnoDB", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/001_ori");
//        this.changeChrPos();
//        this.mergebySub();
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/001_total", "100", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/002_abd", "100", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/003_ab", "100", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/004_d", "100", "1");
        //10 bins
//        this.mkBarplotofDAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/003_daf/001_mkBarplotofDAF", "10", "1");
//        this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/004_maf/001_mkBarplotofMaf","10","0.5");
        //20 bins
//        this.mkBarplotofDAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/005_daf_20bins", "20", "1");
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
            
            
            //开始计算daf
            double[] daf1 = new double[bins];
            TDoubleArrayList dafList1 = new TDoubleArrayList();
            RowTable<String> t = new RowTable<>(infileS);
            int count = t.getRowNumber();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 2).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 2); //DAF值所在的那一列
                dafList1.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf1[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf1.length; i++) {
                daf1[i] = daf1[i] / dafList1.size();
            }
            System.out.println(dafList1.size()  + "  size");
            
            
            //开始计算daf_ABD
            double[] daf2 = new double[bins];
            TDoubleArrayList dafList2 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 3).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 3); //DAF值所在的那一列
                dafList2.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf2[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf2.length; i++) {
                daf2[i] = daf2[i] / dafList2.size();
            }
            System.out.println(dafList2.size() + "  size");
            
            //开始计算daf_AB
            double[] daf3 = new double[bins];
            TDoubleArrayList dafList3 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 4).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 4); //DAF值所在的那一列
                dafList3.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf3[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf3.length; i++) {
                daf3[i] = daf3[i] / dafList3.size();
            }
            System.out.println(dafList3.size() + "  size");
            
            
            //开始写出文件
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Maf\tDensity_Total\tDensity_ABD\tDensity_AB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.3f", (double)bound[i] + (double)(length / bins) / (double)2)).append("\t").append(String.format("%.4f", daf1[i])).append("\t").append(String.format("%.4f", daf2[i])).append("\t").append(String.format("%.4f", daf3[i]));
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
     *
     * @param infileDirS
     * @param outfileDirS
     * @param binNum the number of bins that would be divided
     */
    public void mkBarplotofDAF(String infileDirS, String outfileDirS, String binNum, String max) {
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
            
            
            //开始计算daf
            double[] daf1 = new double[bins];
            TDoubleArrayList dafList1 = new TDoubleArrayList();
            RowTable<String> t = new RowTable<>(infileS);
            int count = t.getRowNumber();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 5).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 5); //DAF值所在的那一列
                dafList1.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf1[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf1.length; i++) {
                daf1[i] = daf1[i] / dafList1.size();
            }
            System.out.println(dafList1.size()  + "  size");
            
            
            //开始计算daf_ABD
            double[] daf2 = new double[bins];
            TDoubleArrayList dafList2 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 6).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 6); //DAF值所在的那一列
                dafList2.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf2[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf2.length; i++) {
                daf2[i] = daf2[i] / dafList2.size();
            }
            System.out.println(dafList2.size() + "  size");
            
            //开始计算daf_AB
            double[] daf3 = new double[bins];
            TDoubleArrayList dafList3 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 7).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 7); //DAF值所在的那一列
                dafList3.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf3[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf3.length; i++) {
                daf3[i] = daf3[i] / dafList3.size();
            }
            System.out.println(dafList3.size() + "  size");
            
            
            //开始写出文件
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Daf\tDensity_Total\tDensity_ABD\tDensity_AB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.3f", (double)bound[i] + (double)(length / bins) / (double)2)).append("\t").append(String.format("%.4f", daf1[i])).append("\t").append(String.format("%.4f", daf2[i])).append("\t").append(String.format("%.4f", daf3[i]));
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
    
    
    public void mergebySub(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub";
        File[] fs = new File(infileDirS).listFiles();
        for(int i = 0; i < fs.length; i++){
            if(fs[i].isHidden())
                fs[i].delete();
        }
        fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }
        for(int i = 0; i < fs.length; i++){
            new CountSites().mergeTxtbysuffix(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath() + "_A.txt.gz", "A.");
            new CountSites().mergeTxtbysuffix(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath()+ "_B.txt.gz", "B.");
            new CountSites().mergeTxtbysuffix(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath()+ "_D.txt.gz", "D.");
        }
        
    }

    public void changeChrPos() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/001_ori";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";
        File[] fs = new File(infileDirS).listFiles();
        for(int i = 0; i < fs.length; i++){
            if(fs[i].isHidden())
                fs[i].delete();
        }
        fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }
        
        for(int i = 0; i < fs.length; i++){
            new CountSites().mergefileandChangeChrPos_chr1and2(fs[i].getAbsolutePath(),new File(outfileDirS, fs[i].getName()).getAbsolutePath());
        }
    }

    /**
     * 将SNP按照同义非同义突变进行分类，并画出DAF分布图
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void classifySNP_byPop(String infileDirS, String outfileDirS) {

        String[] snpClass = {"delSNP", "nonsyTolerantSNP", "synSNP"};
        String[] out = new String[snpClass.length];
        for (int i = 0; i < snpClass.length; i++) {
            new File(outfileDirS, snpClass[i]).mkdirs();
            out[i] = new File(outfileDirS, snpClass[i]).getAbsolutePath();
        }

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            if (fs[i].getName().endsWith(".xlsx")) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String MAF_ABD = null;
                String MAF_AB = null;
                String chrS = f.getName().substring(3, 6);
                boolean ifd = false;
                //根据染色体号进行AB还是D的判断
                String[] db = {"5", "6", "11", "12", "17", "18", "23", "24", "29", "30", "35", "36", "41", "42"};
                Arrays.sort(db);
                if (Arrays.binarySearch(db, chrS) > -1) { //说明是属于D的
                    ifd = true;
                }

                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }

                BufferedWriter[] bw = new BufferedWriter[snpClass.length];
                for (int i = 0; i < bw.length; i++) {
                    outfileS = new File(out[i], "chr" + chrS + "." + snpClass[i] + ".txt.gz").getAbsolutePath();
                    bw[i] = IOUtils.getTextGzipWriter(outfileS);

                    if (ifd == false) {
                        bw[i].write("Chr\tPos\tMaf\tMAF_ABD\tMAF_AB\tDaf\tDaf_ABD\tDaf_AB\tTrans");
                        bw[i].newLine();
                    } else if (ifd == true) {
                        bw[i].write("Chr\tPos\tMaf\tMAF_ABD\tMAF_D\tDaf\tDaf_ABD\tDaf_D\tTrans");
                        bw[i].newLine();
                    }
                }

                String temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
//0Chr	1Pos	2Ref	3Alt	4Major	5Minor	6Maf	7AAF_ABD	8AAF_AB	9Ancestral	10Daf	11Daf_ABD	12Daf_AB	13Variant_type	14SIFT_score	15Transcript
                    List<String> l = PStringUtils.fastSplit(temp);
                    String chr = l.get(0);
                    String pos = l.get(1);
                    String major = l.get(4);
                    String minor = l.get(5);
                    String maf = l.get(6);
                    String AAF_ABD = l.get(7); //注意AAF_ABD中含有0 和1
                    String AAF_AB = l.get(8);
                    String anc = l.get(9);
                    String daf = l.get(10);
                    String Daf_ABD = l.get(11); //注意Daf_ABD中 不含有0和1
                    String Daf_AB = l.get(12);
                    String snpType = l.get(13);
                    String scoreS = l.get(14);
                    String trans = l.get(15);

                    if (AAF_ABD.equals("0.0000") || AAF_ABD.equals("1.0000")) {
                        MAF_ABD = "NA";
                    }
                    if (AAF_AB.equals("0.0000") || AAF_AB.equals("1.0000")) {
                        MAF_AB = "NA";
                    }
                    if (!AAF_ABD.equals("0.0000") && (!AAF_ABD.equals("1.0000"))) { //AAF有分离
                        MAF_ABD = String.valueOf(Math.min(Double.parseDouble(AAF_ABD), 1 - Double.parseDouble(AAF_ABD)));
                    }
                    if (!AAF_AB.equals("0.0000") && (!AAF_AB.equals("1.0000"))) { //AAF有分离
                        MAF_AB = String.valueOf(Math.min(Double.parseDouble(AAF_AB), 1 - Double.parseDouble(AAF_AB)));
                    }

                    //先过滤没有type类型的位点，只保留有类型的位点
                    if (snpType.equals("NA")) {
                        continue;
                    }
                    if (snpType.equals("NONSYNONYMOUS")) { //在类型下进行sift值的判断，
                        if (scoreS.equals("NA")) {
                            continue;
                        }
                        double score = Double.parseDouble(l.get(14));
                        if (score < 0.05) {//说明是有害突变
                            bw[0].write(chr + "\t" + pos + "\t" + maf + "\t" + MAF_ABD + "\t" + MAF_AB + "\t" + daf + "\t" + Daf_ABD + "\t" + Daf_AB + "\t" + trans);
                            bw[0].newLine();
                        } else {//说明是可忍受突变
                            bw[1].write(chr + "\t" + pos + "\t" + maf + "\t" + MAF_ABD + "\t" + MAF_AB + "\t" + daf + "\t" + Daf_ABD + "\t" + Daf_AB + "\t" + trans);
                            bw[1].newLine();
                        }
                    }
                    if (snpType.equals("SYNONYMOUS")) {
                        bw[2].write(chr + "\t" + pos + "\t" + maf + "\t" + MAF_ABD + "\t" + MAF_AB + "\t" + daf + "\t" + Daf_ABD + "\t" + Daf_AB + "\t" + trans);
                        bw[2].newLine();
                    }
                }
                for (int i = 0; i < snpClass.length; i++) {
                    bw[i].flush();
                    bw[i].close();
                }
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileDirS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getCDSannotation(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "chr");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_CDSregion.txt.gz").getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }

            try {
                String header = br.readLine();
                bw.write(header);
                bw.newLine();
                String temp = null;
                List<String> l = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String trans = l.get(15);
                    if (trans.equals("NA")) {
                        continue;
                        /*如果siftascore的值为NA，则无法判断其为有害或是中性突变。我们要筛选即有sift变异类型又有sift值的sites*/
                    }
                    cnt++;
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\t" + String.valueOf(cnt) + "\t" + "trans sites is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * Goal:将ancestral allele添加到数据库中，并计算Daf,Daf_ABD Daf_AB Daf_D
     *
     * @param dbfileS
     * @param ancS
     * @param outfileS
     */
    public void mkSNPsummary_step2(String dbfileS, String ancS, String outfileS) {
        boolean ifd = false;
        double daf = Double.NaN;
        double daf_ABD = Double.NaN;
        double daf_AB = Double.NaN;
        double daf_D = Double.NaN;
        int cntAncNum = 0;
        File f = new File(ancS); //根据ancestral allele 文件，得到染色体号
        int chr = Integer.parseInt(f.getName().substring(3, 6));
        //根据染色体号进行AB还是D的判断
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, chr) > -1) { //说明是属于D的
            ifd = true;
        }
        TIntArrayList snpPosList = new TIntArrayList();
        HashMap<Integer, String> hm = new HashMap<>();

        /*==================================== 建立ancestral allele HashMap =============================================*/
        try { // chr001.wheat.ancestralAllele.txt  chr001_vmap2.1_AnnoDB.txt.gz
            BufferedReader br = null;
            if (f.getName().endsWith(".txt")) {
                br = IOUtils.getTextReader(ancS);
            } else if (f.getName().endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(ancS);
            }
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                cntAncNum++;
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String anc = l.get(3); //
                snpPosList.add(pos);
                hm.put(pos, anc);
            }
            br.close();
            System.out.println(f.getName() + "\tis completed on posList DB with ancestral allele number " + cntAncNum);
        } catch (Exception e) {
            e.printStackTrace();
        }

        int[] snpPos = snpPosList.toArray(new int[snpPosList.size()]);
        Arrays.sort(snpPos);

        try {
            String chrS = PStringUtils.getNDigitNumber(3, chr);
            BufferedReader br = null;
            if (dbfileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(dbfileS);
            } else if (dbfileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(dbfileS);
            }
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }

            String temp = br.readLine(); //read header

            if (ifd == false) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_AB");
                bw.newLine();
            } else if (ifd == true) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_D");
                bw.newLine();
            }

            int cntAnc = 0;
            int cntAncNotMajororMinor = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String major = l.get(4);
                String minor = l.get(5);
                double maf = Double.parseDouble(l.get(6));
                double AAF_ABD = Double.parseDouble(l.get(7));
                double AAF_AB = Double.parseDouble(l.get(8));
                int index = Arrays.binarySearch(snpPos, pos);
                StringBuilder sb = new StringBuilder();
                if (index > -1) { //表明含有anc
                    String ancAllele = hm.get(pos);
                    //如果ancestral allele存在,且等于major，则derived allele等于minor, daf 就等于maf
                    //如果ancestral allele存在,且等于minor，则derived allele等于major, daf 就等于 1-daf1
                    if (ancAllele.equals(minor)) {
                        cntAnc++;
                        daf = 1 - maf;
                        if (AAF_ABD > 0.5) {
                            daf_ABD = AAF_ABD;
                        } else if (AAF_ABD < 0.5) {
                            daf_ABD = 1 - AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = 1 - AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (ancAllele.equals(major)) {
                        cntAnc++;
                        daf = maf;
                        if (AAF_ABD > 0.5) {
                            daf_ABD = 1 - AAF_ABD;
                        } else if (AAF_ABD < 0.5) {
                            daf_ABD = AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = 1 - AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (!ancAllele.equals(minor) && (!ancAllele.equals(major))) {
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                        cntAncNotMajororMinor++;
                    }

                } else { //表明不含anc
                    sb.append(temp).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            double ratio = (double) cntAncNotMajororMinor / (cntAncNotMajororMinor + cntAnc);
            bw.flush();
            bw.close();
            br.close();
            System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
            System.out.println(new File(dbfileS).getName() + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor. The ratio is " + String.format("%.4f", ratio));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 目的：1.将vmap2的chr pos 提取出来，建立数据库；
     *
     * @param infileDirS
     * @param infile2DirS
     * @param outfileDirS
     */
    public void mkSNPsummary_step1(String infileS, String outfileS) {
        //Chr	Pos	Ref	Alt	Major	Minor	Maf	AAF_ABD	AAF_AB
        boolean ifd = false;
        File f = new File(infileS);
        int CHR = Integer.parseInt(f.getName().substring(3, 6));
        //根据染色体号进行AB还是D的判断
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, CHR) > -1) { //说明是属于D的
            ifd = true;
        }

        try {
            BufferedReader br = null;
            BufferedWriter bw = null; // IOUtils.getTextGzipWriter(outfileS);
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }

            if (ifd == false) {
                bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB");
                bw.newLine();
            } else if (ifd == true) {
                bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D");
                bw.newLine();
            }

            String temp = null;
            String te[] = null;
            String major = null;
            String minor = null;
            int biallelicNum = 0;

            String AAF_ABD = null;
            String AAF_AB = null;
            String AAF_D = null;

            while ((temp = br.readLine()) != null) {
                int genoNum = 0;
//                double homNum = 0;
//                double hetNum = 0;
//                double hetRate = 0;
//                double missNum = 0;
//                double missRate = 0;

                double refAlleleGametes = 0;
                double altAlleleGametes = 0;
                double refAF = 0;
                double altAF = 0;
                double maf = 0;
                if (temp.startsWith("#")) {
                    //bw.write(temp);
                    //bw.newLine();
                } else {
                    te = temp.split("\t");
                    String chr = PStringUtils.fastSplit(temp).get(0);
                    String pos = PStringUtils.fastSplit(temp).get(1);
                    String ref = PStringUtils.fastSplit(temp).get(3);
                    String alt = PStringUtils.fastSplit(temp).get(4);

                    if (te[4].length() == 1) { //不含有逗号的情况，即只有一个alt。又开始分，是D 是I 是ATGC 3种情况
                        if (alt.contains("D") || alt.contains("I")) {
                            continue; //只有一个alt且不是indel
                        }
                        AAF_ABD = te[7].split(";")[7].split("=")[1];
                        AAF_AB = te[7].split(";")[8].split("=")[1];
                        AAF_D = te[7].split(";")[8].split("=")[1];

                        biallelicNum++;
                        for (int i = 9; i < te.length; i++) {
                            if (te[i].startsWith(".")) {
//                                missNum++;
                            }
                            if (!te[i].startsWith(".")) {
                                genoNum++; //have the genotype
                                if (te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
//                                    hetNum++; //the number of heterozygous
                                    refAlleleGametes++;
                                    altAlleleGametes++;
                                }
                                if (te[i].startsWith("0/0")) {
//                                    homNum++; //the number of heterozygous
                                    refAlleleGametes++;
                                    refAlleleGametes++;
                                }
                                if (te[i].startsWith("1/1")) {
//                                    homNum++;
                                    altAlleleGametes++;
                                    altAlleleGametes++;
                                }
                            }
                        }
//                        hetRate = hetNum / genoNum;
//                        missRate = missNum / (missNum + genoNum);
                        refAF = refAlleleGametes / (refAlleleGametes + altAlleleGametes);
                        altAF = altAlleleGametes / (refAlleleGametes + altAlleleGametes);;
                        if (refAF > altAF) {
                            major = ref;
                            minor = alt;
                            maf = altAF;
                        } else {
                            maf = refAF;
                            major = alt;
                            minor = ref;
                        }
                        StringBuilder sb = new StringBuilder();
                        //bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB");
                        sb.append(chr).append("\t").append(pos).append("\t").append(ref).append("\t").append(alt).append("\t").
                                append(major).append("\t").append(minor).append("\t").append(String.format("%.4f", maf)).append("\t").
                                append(AAF_ABD).append("\t").append(AAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();

                    }
                }
            }
            System.out.println(infileS + " is completed at " + outfileS);
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
     * 对文件进行分bin，画分布图
     */
    private void mkBarplotOfSNPs() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/delSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/delSNPCount.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/delSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/delSNP/dafSFS.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/nonsyTolerantSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/nonsyTolerantSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/nonsyTolerantSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/nonsyTolerantSNP/dafSFS.txt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/synSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/synSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/synSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/synSNP/dafSFS.txt";
        /**
         * ******************** D subgenome ****************************
         */
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/delSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/delSNPCount.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/delSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/delSNP/dafSFS.txt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/nonsyTolerantSNP/";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/nonsyTolerantSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/nonsyTolerantSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/nonsyTolerantSNP/dafSFS.txt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/synSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/synSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/synSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/synSNP/dafSFS.txt";
        
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub";
        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/synSNP.txt";
        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/mafSFS.txt";
        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/dafSFS.txt";
        
        

        //int sampleSize = 10000;
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                //System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);//将文件删除后，重新将文件列表打印出来，此时，fs不包含隐藏文件。
        /*建立一个边界数组bound，大小是100；bound[i] = 1/100*i;即，均分为100等分！
        建立一个二维数组mafFrequency 和 dafFrequency， 长度为class的种类长，100宽；
        建立一个count 和 dafCount 数组，长度为class分类长；
        建立一个 dafList1 和 dafList 数组，长度为class分类长；
        
        进入for循环，对class文件一一遍历，以读表格的形式进入文件
        Chr	Pos	MinorAllele	MAF	DerivedAllele	DAF
        1	92716	A	0.0077619664	NA	NA
        1	93774	G	8.130081E-4	A	0.9991869919
        1	122123	C	0.0012254902	C	0.0012254902
        做以下几件事情：1，数行数，看每个分类的个数；每个位点一定有maf的值，但不一定有daf的值，因为DA allele不一定存在，若不存在，则DAF值为空。
        2，将maf的值加入mafList;在bound数组里搜索，如果未搜到，则index = - index -2; 否则mafFrequency[i][index]++;
        如果DAF不以N开头，则将daf的值加入 dafList;在bound数组里搜索，如果未搜到，则index = - index -2; 否则dafFrequency[i][index]++; dafCount数组加一
        在此循环内，进入进入另一个for循环j，做统计：
        计算出maf 和daf 在 index 1-100 bin范围内,各个位点所占的百分比。即 index=1, 所有位点count= rownumber， frenquency = index /count；
         */
        int size = 100; //把daf值分成100份，每份有1/100=0.001长度
        double[] bound = new double[size];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double) 1 / size * i;
        }
        double[][] mafFrequency = new double[fs.length][size]; //变异类型的个数
        double[][] dafFrequency = new double[fs.length][size];
        int[] count = new int[fs.length]; //每个变异类型的个数
        int[] dafCount = new int[fs.length];
        TDoubleArrayList[] mafList = new TDoubleArrayList[fs.length]; //每种变异类型的值的集合
        TDoubleArrayList[] dafList = new TDoubleArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            mafList[i] = new TDoubleArrayList();
            dafList[i] = new TDoubleArrayList();
            String infileS = fs[i].getAbsolutePath();
            RowTable<String> t = new RowTable<>(infileS);
            count[i] = t.getRowNumber();
            for (int j = 0; j < t.getRowNumber(); j++) {
                double value = t.getCellAsDouble(j, 2);
                mafList[i].add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                mafFrequency[i][index]++; //值落入第i种变异的第index个区间的个数
                if (!t.getCell(j, 3).startsWith("N")) {
                    value = t.getCellAsDouble(j, 3);
                    dafList[i].add(value);
                    index = Arrays.binarySearch(bound, value);
                    if (index < 0) {
                        index = -index - 2;
                    }
                    dafFrequency[i][index]++;
                    dafCount[i]++;
                }
            }
            for (int j = 0; j < mafFrequency[i].length; j++) {
                mafFrequency[i][j] = mafFrequency[i][j] / count[i];
                dafFrequency[i][j] = dafFrequency[i][j] / dafCount[i]; //因为daf值不是每个都有，有些pos是NA值，所以需要重新计算。
            }
        }
        /*打表格，输入表头，去掉最后一个\t键；表头为4个文件的文件名；
        第二行输入每个分类的conut数；
         */
        try {
            BufferedWriter bw = IOUtils.getTextWriter(countFileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fs.length; i++) {
                sb.append(fs[i].getName().replaceFirst(".changeChrPos.txt", "")).append("\t");
            }
            sb.deleteCharAt(sb.length() - 1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length - 1; i++) {
                bw.write(String.valueOf(count[i]) + "\t");
            }
            bw.write(String.valueOf(count[count.length - 1]));
            bw.newLine();
            bw.flush();
            bw.close();

            bw = IOUtils.getTextWriter(mafDistrubutionFileS); //开始写新的文件MAF
            bw.write("MAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t" + fs[i].getName().replaceFirst(".changeChrPos.txt", "")); //循环写表头
            }
            bw.newLine();
            /*double[][] mafFrequency = new double[fs.length][size] 文件长度为4，size为100*/
            for (int i = 0; i < mafFrequency[0].length; i++) { //i小于第一个文件的长度100，
                bw.write(String.format("%.2f", bound[i]));
                for (int j = 0; j < mafFrequency.length; j++) { //将1-100的频率写出来
                    bw.write("\t" + mafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();

            bw = IOUtils.getTextWriter(dafDistrubutionFileS);//开始写新的文件DAF
            bw.write("DAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t" + fs[i].getName().replaceFirst(".changeChrPos.txt", ""));
            }
            bw.newLine();
            for (int i = 0; i < dafFrequency[0].length; i++) {
                bw.write(String.format("%.2f", bound[i]));
                for (int j = 0; j < dafFrequency.length; j++) {
                    bw.write("\t" + dafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 将SNP按照同义非同义突变进行分类，并画出DAF分布图
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void classifySNPs(String infileDirS, String outfileDirS) {
        String[] snpClass = {"delSNP", "nonsyTolerantSNP", "synSNP"};
        String[] out = new String[snpClass.length];
        for (int i = 0; i < snpClass.length; i++) {
            new File(outfileDirS, snpClass[i]).mkdirs();
            out[i] = new File(outfileDirS, snpClass[i]).getAbsolutePath();
        }

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.stream().forEach(f -> {
            try {
                String chrS = f.getName().substring(3, 6);
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }

                BufferedWriter[] bw = new BufferedWriter[snpClass.length];
                for (int i = 0; i < bw.length; i++) {
                    outfileS = new File(out[i], "chr" + chrS + "." + snpClass[i] + ".txt.gz").getAbsolutePath();
                    bw[i] = IOUtils.getTextGzipWriter(outfileS);
                    bw[i].write("Chr\tPos\tMaf\tDaf\tTrans");
                    bw[i].newLine();
                }
                String temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    String chr = l.get(0);
                    String pos = l.get(1);
                    String major = l.get(4);
                    String minor = l.get(5);
                    double maf = Double.parseDouble(l.get(6));
                    String snpType = l.get(9);
                    String scoreS = l.get(10);
                    String trans = l.get(11);
                    String anc = l.get(12);
                    String daf = l.get(13);
                    //先过滤没有type类型的位点，只保留有类型的位点
                    if (snpType.equals("NA")) {
                        continue;
                    }
                    if (snpType.equals("NONSYNONYMOUS")) { //在类型下进行sift值的判断，
                        if (scoreS.equals("NA")) {
                            continue;
                        }
                        double score = Double.parseDouble(l.get(10));
                        if (score < 0.05) {//说明是有害突变
                            bw[0].write(chr + "\t" + pos + "\t" + maf + "\t" + daf + "\t" + trans);
                            bw[0].newLine();
                        } else {//说明是可忍受突变
                            bw[1].write(chr + "\t" + pos + "\t" + maf + "\t" + daf + "\t" + trans);
                            bw[1].newLine();
                        }
                    }
                    if (snpType.equals("SYNONYMOUS")) {
                        bw[2].write(chr + "\t" + pos + "\t" + maf + "\t" + daf + "\t" + trans);
                        bw[2].newLine();
                    }
                }
                for (int i = 0; i < snpClass.length; i++) {
                    bw[i].flush();
                    bw[i].close();
                }
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileDirS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 解析玉米的基因总结分析
     */
    private void summarizeTranscript2() {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
//        String infileDirS = "/data1/home/aoyue/maizeLoad/001_variantSummary/001_hmp321Info_filter";
//        String geneFeatureFileS = "/data1/publicData/maize/gene/Zea_mays.AGPv4.38.pgf";
//        String outfileS = "/data1/home/aoyue/maizeLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
        double gerpCut = 0;
        File[] fs = new File(infileDirS).listFiles();
        int chrNum = fs.length;
        HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
        int[][] snpPos = new int[chrNum][];
        byte[][] snps = new byte[chrNum][];
        byte[][] snpAnc = new byte[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            posGeneMap[i] = new HashMap();
        }
        //下面这一段将posGeneMap建立完整，使每个位点对应哪些基因名字，都装进这个map里 
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        HashMap<String, Integer> geneCDSLengthMap = new HashMap();
        /*将所有基因的名字进行for循环输入到数组genes中，对应于每一个基因，我们通过getTranscriptName得到转录本的名字，通过getCDSList方法得到编码序列的起始位点*/
        List<String> genesList = new ArrayList<>();
        //String[] genes = new String[gf.getGeneNumber()]; //这个是原来的代码，后续需要修改，因为gf文件中包含我们不需要的染色体上的基因
        int cntchr11and12 = 0;
        int cntchr1to10 = 0;

        //*********************************** START1 ***********************************//
        //该段代码的作用是，通过读取每个基因，得到最长转录本的名字，计算该转录本的长度。
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chrIndex = gf.getGeneChromosome(i) - 1;
            /*这个地方是先过滤数据，将定位在11号12号染色体上的基因过滤掉，并且跳出循环*/
            if (chrIndex > 9) {
                cntchr11and12++;
                continue;
            }
            cntchr1to10++; //能够得到1-10号染色体的基因数目
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String geneName = gf.getTranscriptName(i, longTransIndex); //得到最长的转录本的名字
            //genes[i] = geneName;
            genesList.add(geneName);
            List<Range> cdsList = gf.getCDSList(i, longTransIndex);
            /*得到基因的最长转录本的CDSList*/
            int cnt = 0;

            /*对于每一个基因的编码序列，还有很多个cds片段，即cdsList；我们对cdsList进行for循环，得到每个cds的起始和终止位置，从而计算出总长*/
            for (int j = 0; j < cdsList.size(); j++) {
                int rStart = cdsList.get(j).start;
                int rEnd = cdsList.get(j).end;
                for (int k = rStart; k < rEnd; k++) {
                    /*posGeneMap是一个HashMap数组，一条染色体对应一个String类型的ArrayList；
                    故得到该位点的所属基因名字列表，如果该位点不含基因名，就将 genename赋值给该位点，完善posGeneMap
                    否则，如果该位点含有其他基因的基因名字，依旧把genename赋值给该位点*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(k); //建立map的关系，那个位点对应哪个list HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum]; 
                    if (geneNameList == null) {
                        geneNameList = new ArrayList();
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                    } else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                        /*最终将posGeneMap绘图完成*/
                    }
                    cnt++;
                    /*每一个CDS位点相加，最终得到这个cds的长度。*/
                }

                // 最终cnt是一个基因的所有cdslist中，每个cds的每个位点包含的基因数目的总和
            } //该循环是一个基因的所有cds循环
            geneCDSLengthMap.put(geneName, cnt); //
        }
        //*********************************** END1 ***********************************//

        System.out.println(cntchr11and12 + "genes are not used");
        System.out.println(cntchr1to10 + "genes are used");

        String[] genes = genesList.toArray(new String[genesList.size()]);
        Arrays.sort(genes); //genes指所有基因对应最长转录本的名字的组合，是一个数组。
        int[] snpCount = new int[genes.length]; //这一步也很重要，就是想知道每个基因内部有多少个变异位点。
        List<File> hmpList = Arrays.asList(fs); //fs指的是根据VCF位点建立的注释数据库，该数据库包含很多计算的信息，有SIFT值GERP值。
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", "")) - 1;
            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            TByteArrayList snpAncList = new TByteArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                /*这一段主要是将有基因的位点列出来，及将转录组的位点列出来*/
                while ((temp = br.readLine()) != null) {
                    cnt++; //检测程序运行的情况
                    if (cnt % 1000000 == 0) {
                        System.out.println("Hmp\tchr" + String.valueOf(chrIndex + 1) + "\t" + String.valueOf(cnt) + " ###hmpInfo Process");
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).contains("<") || l.get(3).contains(",")) {
                        continue; //过滤含有2个alt和含有indel的位点
                    }
                    int pos = Integer.valueOf(l.get(1));
                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) {
                        continue; //说明该变异位点不在基因区
                    }
                    for (int i = 0; i < geneNameList.size(); i++) {
                        int index = Arrays.binarySearch(genes, geneNameList.get(i));//在基因库的第i个位置，该基因数加一
                        snpCount[index]++; //第i个位置的基因含有的snp数目；
                    }
                    /*将该位置加入snpPosList，表明该位置是有基因的，没有基因的位置就不必加到snpPosList中去*/
                    snpPosList.add(pos);
                    snpList.add(l.get(3).getBytes()[0]); //返回该位点的alt碱基的二进制码
                    snpAncList.add(l.get(4).getBytes()[0]);  //返回该位点的ancestral碱基的二进制码
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            snpPos[chrIndex] = snpPosList.toArray(); //第n条染色体的在基因区间的snp所对应的pos的集合；
            snps[chrIndex] = snpList.toArray(); //第n条染色体的在区间内的snp所对应的alt的集合
            snpAnc[chrIndex] = snpAncList.toArray(); //第n条染色体的在区间内的snp所对应的ancestral allele的集合
        });

        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length]; //非同义突变的数目
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];  //非同义突变，但是sift值是NA的数目

        int[] b73SynCount = new int[genes.length];
        int[] b73NonCount = new int[genes.length];
        int[] b73DelCount = new int[genes.length];
        int[] b73DelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];

        TIntArrayList[] delPosList = new TIntArrayList[chrNum]; //有害变异的位点集合
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", "")) - 1;
            delPosList[chrIndex] = new TIntArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt % 1000000 == 0) {
                        System.out.println("Sift\tchr" + String.valueOf(chrIndex + 1) + "\t" + String.valueOf(cnt) + " ### SIFT Process");
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(16).startsWith("NA")) {
                        continue; //Variant_type	SIFT_score	Transcript  16 17 18
                    }
                    if (l.get(18).startsWith("NA")) {
                        continue; //没有变异类型和转录本的位点，都过滤掉。
                    }
                    String gene = l.get(18);
                    int geneIndex = Arrays.binarySearch(genes, gene); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) {
                        continue;
                    }
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos); // 
                    if (index < 0) {
                        continue;
                    }
                    if (snps[chrIndex][index] != l.get(3).getBytes()[0]) {
                        continue; //再次验证，如果该位点的alt和数据库中的alt不一致，则过滤掉。
                    }
                    byte ref = l.get(2).getBytes()[0];
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (snpAnc[chrIndex][index] == snps[chrIndex][index]) { //如果ancestral allele 等于alt的话，derived allele就等于1
                        derivedState = 1; //mean b73 carries derived allele
                    } else if (snpAnc[chrIndex][index] == ref) { //如果ancestral allele 等于ref的话，derived allele就等于0
                        derivedState = 0;
                    }

                    if (derivedState == -1) {
                        noAncCount[geneIndex]++; //如果ancestral allele 不存在的话，derived allele就等于-1
                    }
                    String type = null;
                    if (l.get(16).equals("NA")) {

                    } else {
                        if (l.get(16).equals("SYNONYMOUS")) { //如果type等于syn，那么 该位点所属的基因的syn属性就加一
                            type = "Syn";
                            synCount[geneIndex]++;
                            if (derivedState == 1) {  //如果derived allele就等于1，说明ancestral allele 等于alt， derived allele 等于ref； 如何计算daf,判断da是major还是minor，如果da是major那么daf=1-daf1，如何判断major allele和minor allele？ 如果ref allele frequency > alt allele frequence,那么major是ref; 反之亦然；
                                b73SynCount[geneIndex]++; //如果参考基因组是 derived allele,那么就加一，为什么？
                            }
                        } else {
                            type = "Non";
                            nonCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73NonCount[geneIndex]++;
                            }
                            if (l.get(17).startsWith("NA")) { //index 17列是sift的值，NON-SYNONYMOUS 存在的情况下，sift可能有也可能没有。
                                naCount[geneIndex]++; //是nonsynonymous类型但是没有sift值的个数
                            } else {
                                if (Double.valueOf(l.get(17)) < 0.05) {
                                    delCount[geneIndex]++;
                                    delPosList[chrIndex].add(pos);
                                    if (derivedState == 1) {
                                        b73DelCount[geneIndex]++;
                                    }
                                }
                            }
                        }
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        /**
         *
         */
        int[][] delPos = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delPos[i] = delPosList[i].toArray();
            Arrays.sort(delPos[i]);
        }
        int[] gerpAlignCount = new int[genes.length];
        int[] snpGerpAlignCount = new int[genes.length];
        double[] gerpTree = new double[genes.length];
        double[] gerpScore = new double[genes.length];
        double[] snpGerpTree = new double[genes.length];
        double[] snpGerpScore = new double[genes.length];

        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", "")) - 1;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) { //gerp文件没有表头
                    cnt++;
                    if (cnt % 1000000 == 0) {
                        System.out.println("Gerp\tchr" + String.valueOf(chrIndex + 1) + "\t" + String.valueOf(cnt) + " ### Gerp Process");
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(1));
                    if (l.get(14).equals("NA")) {
                        continue;
                    }
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos); //根据pos信息，得到该pos对应的gene name的集合
                    if (geneNameList == null) {
                        continue;
                    }
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        if (gene == null) {
                            continue;
                        }
                        int geneIndex = Arrays.binarySearch(genes, gene);

                        double treeValue = Double.valueOf(l.get(14));
                        double scoreValue = Double.valueOf(l.get(15));
                        if (treeValue == 0) {
                            continue; //过滤枝长是0的数目
                        }
                        gerpAlignCount[geneIndex]++; //如果枝长不是0，说明该位点存在保守不保守
                        gerpTree[geneIndex] += treeValue;
                        gerpScore[geneIndex] += scoreValue; //第i个基因的gerpscore的总和是多少
                        int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                        if (index < 0) {
                            continue;
                        }

                        byte derivedState = 0; //mean ancestral allele is not defined or ancestral allele is alt
                        if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                            derivedState = 1; //mean b73 carries derived allele
                        }
                        snpGerpAlignCount[geneIndex]++;
                        snpGerpTree[geneIndex] += treeValue;
                        snpGerpScore[geneIndex] += scoreValue;
                        index = Arrays.binarySearch(delPos[chrIndex], pos);
                        if (index < 0) {
                            continue;
                        }
                        if (scoreValue <= gerpCut) {
                            continue;
                        }
                        delHGCount[geneIndex]++;
                        if (derivedState == 1) {
                            b73DelHGCount[geneIndex]++;
                        }
                    }

                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        gf.sortGeneByStartPosition();
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpTreeLength\tMeanGerpScore\tSNPMeanGerpTreeLength\tSNPMeanGerpScore";
            header = header + "\tNumAmbigousAnc\tB73NumberOfSyn\tB73PercentageSyn\tB73NumberOfNon\tB73PercentageNon\tB73NumberOfDeleterious\tB73PercentageDeleterious\tB73NumberOfHGDeleterious\tB73PercentageHGDeleterious";
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < genes.length; i++) {
                //if(gf.getGeneChromosome(i) > 10) continue;
                StringBuilder sb = new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append((double) snpCount[i] / cdsLength).append("\t");
                int ifSiftAligned = 1;
                if (naCount[i] == nonCount[i]) {
                    ifSiftAligned = 0; //非同义突变，但是sift值是NA的数目 等于 非同义突变的数目，那么说明在这个基因内部没有sift突变
                }
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append((double) synCount[i] / cdsLength).append("\t");
                sb.append(nonCount[i]).append("\t").append((double) nonCount[i] / cdsLength).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) {
                    ratio = Double.NaN;
                } else {
                    ratio = (double) nonCount[i] / synCount[i];
                }
                sb.append(ratio).append("\t").append(delCount[i]).append("\t").append((double) delCount[i] / cdsLength).append("\t");
                sb.append(delHGCount[i]).append("\t").append((double) delHGCount[i] / cdsLength).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) {
                    ifGerpAligned = 0;
                }
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append((double) gerpAlignCount[i] / cdsLength).append("\t");

                if (gerpAlignCount[i] == 0) {
                    sb.append(Double.NaN).append("\t").append(Double.NaN).append("\t");
                } else {
                    sb.append((double) gerpTree[i] / gerpAlignCount[i]).append("\t").append((double) gerpScore[i] / gerpAlignCount[i]).append("\t");
                }

                if (snpGerpAlignCount[i] == 0) {
                    sb.append(Double.NaN).append("\t").append(Double.NaN);
                } else {
                    sb.append((double) snpGerpTree[i] / snpGerpAlignCount[i]).append("\t").append((double) snpGerpScore[i] / snpGerpAlignCount[i]);
                }

                double cdsL = (double) (snpCount[i] - noAncCount[i]) / snpCount[i] * cdsLength;
                sb.append("\t").append(noAncCount[i]).append("\t").append(b73SynCount[i]).append("\t").append((double) b73SynCount[i] / cdsL).append("\t");
                sb.append(b73NonCount[i]).append("\t").append((double) b73NonCount[i] / cdsL).append("\t");

                sb.append(b73DelCount[i]).append("\t").append((double) b73DelCount[i] / cdsL).append("\t");
                sb.append(b73DelHGCount[i]).append("\t").append((double) b73DelHGCount[i] / cdsL);

                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void scriptAddAncAllele() {
//        for (int i = 1; i < 43; i++) {
//            String CHR = PStringUtils.getNDigitNumber(3, i);
//            System.out.println("java -Xms50g -Xmx100g -jar 017_mkAnnoDB.addAncAllele.single.jar /data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/chr" 
//                    + CHR + ".lineage.maf0.005.bi.AnnoDB.txt.gz "
//                    + "/data4/home/aoyue/vmap2/analysis/ancestralAllele/Chr" + CHR 
//                    + ".ancestralAllele.txt " 
//                    + "/data4/home/aoyue/vmap2/analysis/015_annoDB/002_addAncestralAllele/chr"
//                    + CHR + ".lineage.maf0.005.bi.AnnoDB.addAncAllele.txt.gz"
//                    );
//        }

        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(l, i);
            if (index > -1) { //如果大于-1，则在集合中搜索到，说明是属于D基因组的；
                System.out.println("java -Xms50g -Xmx100g -jar 017_mkAnnoDB.addAncAllele.single.jar /data4/home/aoyue/vmap2/analysis/015_annoDB/003_addSIFT/chr"
                        + chr + ".lineage.maf0.005.bi.AnnoDB.addSIFT.txt.gz "
                        + "/data4/home/aoyue/vmap2/daxing/ancestralAllele/chr" + chr
                        + ".wheat.ancestralAllele.txt "
                        + "/data4/home/aoyue/vmap2/analysis/015_annoDB/005_addAncestralAllele/chr"
                        + chr + ".lineage.maf0.005.bi.AnnoDB.addSIFT.addAnc.txt.gz > log_017/log_" + chr + "_mkAnnoDB.addAncAllele.single.txt"
                );
            }
        }
    }

    /**
     *
     * @param dbfileS
     * @param ancS
     * @param outfileS
     */
    public void addAncAllele_singlethread(String dbfileS, String ancS, String outfileS) {
        double daf = Double.NaN;

        //Step1:建立42个文件的snpPos的二维数组，第一维是染色体号，第二维是每条染色体含有的pos集合；
        //Step2:然后文件并行流读进去，将染色体号提取出来，提取的chr的数字形式-1就是snpPos[][]的第一维；
        //在文件流内建立一个集合posList，将pos信息加入posList中，文件读入完毕后，将该posList转换成数组，放入每一维染色体对应的pos信息。
        //在文件流内建立一个HashMap集合，将pos对应的ancestral allele建立联系，
        //Step3:读入数据库文件，进行一行一行搜索，如果posList中有该位点，说明存在ancestral alle,即在文件最后一列添加，如果没有搜到，则输出NA
        File f = new File(ancS);
        int chr = Integer.parseInt(f.getName().substring(3, 6));
        TIntArrayList snpPosList = new TIntArrayList();
        HashMap<Integer, String> hm = new HashMap<>();

        try { //Chr001.ancestralAllele.txt  chr001.lineage.maf0.005.bi.AnnoDB.txt.gz
            BufferedReader br = null;
            if (f.getName().endsWith(".txt")) {
                br = IOUtils.getTextReader(ancS);
            } else if (f.getName().endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(ancS);
            }
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String anc = l.get(3);
                snpPosList.add(pos);
                hm.put(pos, anc);
            }
            br.close();
            System.out.println(f.getName() + "\tis completed on posList DB ");
        } catch (Exception e) {
            e.printStackTrace();
        }

        int[] snpPos = snpPosList.toArray(new int[snpPosList.size()]);
        Arrays.sort(snpPos);

        try { //chr001.lineage.maf0.005.bi.AnnoDB.txt.gz
            String chrS = PStringUtils.getNDigitNumber(3, chr);
            BufferedReader br = null;
            if (dbfileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(dbfileS);
            } else if (dbfileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(dbfileS);
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //压缩格式的输出
            String temp = br.readLine();
            bw.write(temp + "\tAncestralAllele\tDaf");
            bw.newLine();
            int cntAnc = 0;
            int cntAncNotMajororMinor = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String major = l.get(4);
                String minor = l.get(5);
                double maf = Double.parseDouble(l.get(6));
                int index = Arrays.binarySearch(snpPos, pos);
                StringBuilder sb = new StringBuilder();
                if (index > -1) { //表明含有anc
                    String ancAllele = hm.get(pos);
                    //如果ancestral allele存在,且等于major，则da等于的minor, daf 就等于maf
                    //如果ancestral allele存在,且等于minor，则da等于的major, daf 就等于 1-daf1
                    if (ancAllele.equals(minor)) {
                        cntAnc++;
                        daf = 1 - maf;
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (ancAllele.equals(major)) {
                        cntAnc++;
                        daf = maf;
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (!ancAllele.equals(minor) && !ancAllele.equals(major)) {
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                        cntAncNotMajororMinor++;
                    }

                } else { //表明不含anc
                    sb.append(temp).append("\t").append("NA").append("\t").append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
            System.out.println("chr" + PStringUtils.getNDigitNumber(3, chr) + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor.");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void addAncestralAllele(String infileDirS, String ancestralDirS, String outfileDirS) {
        //String infileDirS = "";
        //String outfileDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/002_addAncestralAllele";

        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }

        //Step1:建立42个文件的snpPos的二维数组，第一维是染色体号，第二维是每条染色体含有的pos集合；
        //Step2:然后文件并行流读进去，将染色体号提取出来，提取的chr的数字形式-1就是snpPos[][]的第一维；
        //在文件流内建立一个集合posList，将pos信息加入posList中，文件读入完毕后，将该posList转换成数组，放入每一维染色体对应的pos信息。
        //在文件流内建立一个HashMap集合，将pos对应的ancestral allele建立联系，
        //Step3:读入数据库文件，进行一行一行搜索，如果posList中有该位点，说明存在ancestral alle,即在文件最后一列添加，如果没有搜到，则输出NA
        File[] fs = new File(ancestralDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(ancestralDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            int chr = Integer.parseInt(f.getName().substring(3, 6));
            TIntArrayList snpPosList = new TIntArrayList();
            HashMap<Integer, String> hm = new HashMap<>();

            try { //Chr001.ancestralAllele.txt  chr001.Alineage.maf0.005.bi.AnnoDB.txt.gz
                String infileS = f.getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String anc = l.get(3);
                    snpPosList.add(pos);
                    hm.put(pos, anc);
                }
                br.close();
                System.out.println(f.getName() + "\tis completed on posList DB ");
            } catch (Exception e) {
                e.printStackTrace();
            }

            int[] snpPos = snpPosList.toArray(new int[snpPosList.size()]);
            Arrays.sort(snpPos);

            try { //chr001.Alineage.maf0.005.bi.AnnoDB.txt.gz
                String chrS = PStringUtils.getNDigitNumber(3, chr);
                String infileS = new File(infileDirS, "chr" + String.valueOf(chrS) + "." + hml.get(chr) + "lineage.maf0.005.bi.AnnoDB.txt.gz").getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, "chr" + String.valueOf(chrS) + "." + hml.get(chr) + "lineage.maf0.005.bi.AnnoDB.addAncAllele.txt.gz").getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, "chr" + String.valueOf(chrS) + "." + hml.get(chr) + "lineage.maf0.005.bi.AnnoDB.addAncAllele.txt.gz").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //压缩格式的输出
                String temp = br.readLine();
                bw.write(temp + "\tAncestralAllele\tDaf");
                bw.newLine();
                int cntAnc = 0;
                int cntAncNotMajororMinor = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String major = l.get(4);
                    String minor = l.get(5);
                    double maf = Double.parseDouble(l.get(6));
                    int index = Arrays.binarySearch(snpPos, pos);
                    if (index > -1) { //表明含有anc
                        String ancAllele = hm.get(pos);
                        double daf = Double.NaN;
                        //如果ancestral allele存在,且等于major，则da等于的minor, daf 就等于maf
                        //如果ancestral allele存在,且等于minor，则da等于的major, daf 就等于 1-daf1
                        if (ancAllele.equals(minor)) {
                            cntAnc++;
                            daf = 1 - maf;
                            StringBuilder sb = new StringBuilder();
                            sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                        if (ancAllele.equals(major)) {
                            cntAnc++;
                            daf = maf;
                            StringBuilder sb = new StringBuilder();
                            sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                        if (!ancAllele.equals(minor) && !ancAllele.equals(major)) {
                            StringBuilder sb = new StringBuilder();
                            sb.append(temp).append("\t").append(ancAllele).append("\t").append("NA");
                            bw.write(sb.toString());
                            bw.newLine();
                            //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                            // cntAncNotMajororMinor

                        }

                    } else if (index < 0) { //表明不含anc
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp).append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
                System.out.println("chr" + PStringUtils.getNDigitNumber(3, chr) + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor.");
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     *
     *
     * @param infileDirS
     */
    public void mkSNPsummary(String infileDirS, String outfileDirS) {

        new File(outfileDirS).mkdirs();
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        /**
         * *************************************************************
         */
        System.out.println("FileName\tbiallelicNum\tBiallelicMafmore0.005Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".AnnoDB.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".AnnoDB.txt.gz")).getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tHetProportion\tMissProportion");
                bw.newLine();
                String temp = null;
                String te[] = null;
                int biallelicNum = 0;
                int biallelicMafmoreNum = 0;

                String major = null;
                String minor = null;
                while ((temp = br.readLine()) != null) {
                    int genoNum = 0;
                    double homNum = 0;
                    double hetNum = 0;
                    double hetRate = 0;
                    double missNum = 0;
                    double missRate = 0;

                    double refAlleleGametes = 0;
                    double altAlleleGametes = 0;
                    double refAF = 0;
                    double altAF = 0;
                    double maf = 0;

                    if (temp.startsWith("#")) {
                        //bw.write(temp);
                        //bw.newLine();
                    } else {
                        te = temp.split("\t");
                        String chr = PStringUtils.fastSplit(temp).get(0);
                        String pos = PStringUtils.fastSplit(temp).get(1);
                        String ref = PStringUtils.fastSplit(temp).get(3);
                        String alt = PStringUtils.fastSplit(temp).get(4);

                        if (te[4].length() == 1) { //不含有逗号的情况，即只有一个alt。又开始分，是D 是I 是ATGC 3种情况
                            if (alt.contains("D") || alt.contains("I")) {
                                continue; //只有一个alt且不是indel
                            }
                            biallelicNum++;
                            for (int i = 9; i < te.length; i++) {
                                if (te[i].startsWith(".")) {
                                    missNum++;
                                }
                                if (!te[i].startsWith(".")) {
                                    genoNum++; //have the genotype
                                    if (te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                        hetNum++; //the number of heterozygous
                                        refAlleleGametes++;
                                        altAlleleGametes++;
                                    }
                                    if (te[i].startsWith("0/0")) {
                                        homNum++; //the number of heterozygous
                                        refAlleleGametes++;
                                        refAlleleGametes++;
                                    }
                                    if (te[i].startsWith("1/1")) {
                                        homNum++;
                                        altAlleleGametes++;
                                        altAlleleGametes++;
                                    }
                                }
                            }
                            hetRate = hetNum / genoNum;
                            missRate = missNum / (missNum + genoNum);
                            refAF = refAlleleGametes / (refAlleleGametes + altAlleleGametes);
                            altAF = altAlleleGametes / (refAlleleGametes + altAlleleGametes);;
                            if (refAF >= altAF) {
                                major = ref;
                                minor = alt;
                                maf = altAF;
                            } else {
                                maf = refAF;
                                major = alt;
                                minor = ref;
                            }

                            if (maf <= 0.005) {
                                continue;
                            }
                            biallelicMafmoreNum++;
                            StringBuilder sb = new StringBuilder();
                            //bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tHetProportion\tMissProportion");
                            sb.append(chr).append("\t").append(pos).append("\t").append(ref).append("\t").append(alt).append("\t").
                                    append(major).append("\t").append(minor).append("\t").append(String.format("%.4f", maf)).append("\t").
                                    append(String.format("%.4f", hetRate)).append("\t").append(String.format("%.4f", missRate));
                            bw.write(sb.toString());
                            bw.newLine();

                        }
                    } //else的终止
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(String.valueOf(f.getName()) + "\t" + String.valueOf(biallelicNum) + "\t" + String.valueOf(biallelicMafmoreNum) + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

}
