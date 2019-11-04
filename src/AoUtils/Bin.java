/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class Bin {

    public Bin() {
        //this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/002_calMAF","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/003_bintable","25","0.5");

//        this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/009_calMAF_newData","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/010_bintable","25","0.5");
//    this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/012_calMAF_bySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/013_bintable_bySub","25","0.5");
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
            if (t.getCell(i, 1).equals("NA")) {
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
        System.out.println(mafList.size()  + "  size");
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
                if (t.getCell(i, 0).equals("NA")) { //MAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 3); //MAF值所在的那一列
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
                System.out.println(f.getName() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * make bin table, 一个文件制作一个Binfile，有起始和终止位置
     *
     * @param infileDirS
     * @param outfileDirS
     * @param binNum
     */
    public void mkBintable(String infileDirS, String outfileDirS, String binNum) {
        int binSize = Integer.parseInt(binNum);
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
            for (int i = 0; i < bound.length; i++) {
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
