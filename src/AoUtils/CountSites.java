/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class CountSites {

    public CountSites() {

    }

    /**
     * 将改变位置的chr5和chr6文本文件合并成1D一个文件， chr11,chr12两个文件合并成2D一个文件。自动识别染色体序号并进行合并。
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mergefile1and2(String infileDirS, String outfileDirS) {

        //建立1-44一一对应chr1A的关系,目的：根据chr1找到chr1A
        String[] chrs = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Mit", "Chl"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44};
        HashMap<Integer, String> hmcntchr = new HashMap<>();

        int cnt = 0;
        for (int i = 0; i < chrs.length; i++) {
            if (cnt == 43) {
                hmcntchr.put(43, "Mit");
            }
            if (cnt == 44) {
                hmcntchr.put(44, "Chl");
            } else {
                hmcntchr.put(cnts[cnt], chrs[i]);
                cnt++;
                hmcntchr.put(cnts[cnt], chrs[i]);
                cnt++;
            }
        }

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);

        try {
            for (int i = 0; i < fs.length; i++) {
                int chr = Integer.valueOf(fs[i].getName().substring(3, 6)); //先对文件的题目进行处理，获取染色体号，进行判断
                for (int j = 1; j < 43; j++) {
                    if (!(j % 2 == 0)) { // j只进行奇数判断，如只进行 chr1 3 5 7 9判断
                        ///******************内部开始写*******************************//
                        if (chr == j) {
                            //读入文件
                            String infileS = fs[i].getAbsolutePath();
                            BufferedReader br = null;
                            if (infileS.endsWith(".txt")) {
                                br = IOUtils.getTextReader(infileS);
                            } else if (infileS.endsWith(".txt.gz")) {
                                br = IOUtils.getTextGzipReader(infileS);
                            }

                            //确定输出文件的路径，并读入header
                            String secondchr = PStringUtils.getNDigitNumber(3, chr + 1);
                            String outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + ".ChrPos.txt.gz").getAbsolutePath();
                            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                            bw.write(br.readLine()); //先读表头
                            bw.newLine();
                            ///开始合并文件1和2

                            String temp = null; //read header
                            while ((temp = br.readLine()) != null) {
                                StringBuilder sb = new StringBuilder();
                                sb.append(temp);
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                            int a =3;
                            //开始读入第2个文件
                            infileS = new File (fs[i].getParent(),fs[i].getName().replaceFirst(PStringUtils.getNDigitNumber(3, chr), secondchr)).getAbsolutePath();
                            if (infileS.endsWith(".txt")) {
                                br = IOUtils.getTextReader(infileS);
                            } else if (infileS.endsWith(".txt.gz")) {
                                br = IOUtils.getTextGzipReader(infileS);
                            }

                            temp = br.readLine(); //read header
                            while ((temp = br.readLine()) != null) {
                                StringBuilder sb = new StringBuilder();
                                sb.append(temp);
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                            br.close();
                            bw.flush();
                            bw.close();

                        }
                        ///********************内部写出结束*****************************//
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 目的：将所有txt文本的chr pos位点信息合并成一个文件。
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergeTxt(String infileDirS, String outfileS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        System.out.println("Chr\tSNP_Num");

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();

            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = br.readLine(); //read header
                //int chrint = Integer.parseInt(fs[i].getName().substring(3, 6));
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
     * Change the chr and pos of vcf file to chr1A 1B pattern.
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void changechrPosOnVCF(String infileDirS, String outfileDirS) {

        //建立1-44一一对应chr1A的关系,目的：根据chr1找到chr1A
        String[] chrs = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Mit", "Chl"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44};
        HashMap<Integer, String> hmcntchr = new HashMap<>();

        int cnt = 0;
        for (int i = 0; i < chrs.length; i++) {
            if (cnt == 43) {
                hmcntchr.put(43, "Mit");
            }
            if (cnt == 44) {
                hmcntchr.put(44, "Chl");
            } else {
                hmcntchr.put(cnts[cnt], chrs[i]);
                cnt++;
                hmcntchr.put(cnts[cnt], chrs[i]);
                cnt++;
            }
        }

        //建立 chr 2 4 6 8 坐标位置修改的hashMap关系，目的：根据chr2,改变chr2的pos信息和chr信息；
        HashMap<Integer, Integer> hmcntsPOS = new HashMap<>();
        int A1 = 471304005;
        int B1 = 438720154;
        int D1 = 452179604;
        int A2 = 462376173;
        int B2 = 453218924;
        int D2 = 462216879;
        int A3 = 454103970;
        int B3 = 448155269;
        int D3 = 476235359;
        int A4 = 452555092;
        int B4 = 451014251;
        int D4 = 451004620;
        int A5 = 453230519;
        int B5 = 451372872;
        int D5 = 451901030;
        int A6 = 452440856;
        int B6 = 452077197;
        int D6 = 450509124;
        int A7 = 450046986;
        int B7 = 453822637;
        int D7 = 453812268;
        int[] chrpos = {A1, B1, D1, A2, B2, D2, A3, B3, D3, A4, B4, D4, A5, B5, D5, A6, B6, D6, A7, B7, D7};
        int tt = 0;
        for (int i = 0; i < cnts.length; i++) {
            if (cnts[i] == 43 || (cnts[i] == 44)) {

            } else {
                if (!(i % 2 == 0)) {
                    hmcntsPOS.put(cnts[i], chrpos[tt]);
                    tt++;
                }
            }
        }

        //多线程处理文件，修改坐标,重新输入文件
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".changeChrPos.vcf.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".changeChrPos.vcf.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int a = 0;
                List<String> l = null;
                int chr = 0;
                int pos = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        l = PStringUtils.fastSplit(temp);
                        StringBuilder sb = new StringBuilder();
                        chr = Integer.valueOf(l.get(0));
                        pos = Integer.valueOf(l.get(1));
                        String Chr = hmcntchr.get(chr);
                        if (chr == 43 || chr == 44) {

                        } else {
                            if (chr % 2 == 0) {
                                pos = pos + hmcntsPOS.get(chr);
                            }
                            sb = new StringBuilder();
                            sb.append(Chr).append("\t").append(String.valueOf(pos));
                            for (int j = 2; j < l.size(); j++) {
                                sb.append("\t");
                                sb.append(l.get(j));
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                            if (a % 1000000 == 0) {
                                System.out.println("Output " + String.valueOf(a) + " SNPs");
                            }
                            a++;
                        }
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(a) + " SNPs output from " + f.getAbsolutePath());
                System.out.println(chr + " was finished to change the CHROM and Pos from    " + f.getAbsolutePath());
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * Change the txt file of chr pos column into 1A mode. Chr Pos 6	291 ---> 1D
     * 452179895
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void changechrPos(String infileDirS, String outfileDirS) {

        //建立1-44一一对应chr1A的关系,目的：根据chr1找到chr1A
        String[] chrs = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Mit", "Chl"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44};
        HashMap<Integer, String> hmcntchr = new HashMap<>();

        int cnt = 0;
        for (int i = 0; i < chrs.length; i++) {
            if (cnt == 43) {
                hmcntchr.put(43, "Mit");
            }
            if (cnt == 44) {
                hmcntchr.put(44, "Chl");
            } else {
                hmcntchr.put(cnts[cnt], chrs[i]);
                cnt++;
                hmcntchr.put(cnts[cnt], chrs[i]);
                cnt++;
            }
        }

        //建立 chr 2 4 6 8 坐标位置修改的hashMap关系，目的：根据chr2,改变chr2的pos信息和chr信息；
        HashMap<Integer, Integer> hmcntsPOS = new HashMap<>();
        int A1 = 471304005;
        int B1 = 438720154;
        int D1 = 452179604;
        int A2 = 462376173;
        int B2 = 453218924;
        int D2 = 462216879;
        int A3 = 454103970;
        int B3 = 448155269;
        int D3 = 476235359;
        int A4 = 452555092;
        int B4 = 451014251;
        int D4 = 451004620;
        int A5 = 453230519;
        int B5 = 451372872;
        int D5 = 451901030;
        int A6 = 452440856;
        int B6 = 452077197;
        int D6 = 450509124;
        int A7 = 450046986;
        int B7 = 453822637;
        int D7 = 453812268;
        int[] chrpos = {A1, B1, D1, A2, B2, D2, A3, B3, D3, A4, B4, D4, A5, B5, D5, A6, B6, D6, A7, B7, D7};
        int tt = 0;
        for (int i = 0; i < cnts.length; i++) {
            if (cnts[i] == 43 || (cnts[i] == 44)) {

            } else {
                if (!(i % 2 == 0)) {
                    hmcntsPOS.put(cnts[i], chrpos[tt]);
                    tt++;
                }
            }
        }

        //多线程处理文件，修改坐标,重新输入文件
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", ".changeChrPos.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt.gz", ".changeChrPos.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = br.readLine();
                bw.write(temp);
                bw.newLine();
                int a = 0;
                List<String> l = null;
                int chr = 0;
                int pos = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    chr = Integer.valueOf(l.get(0));
                    pos = Integer.valueOf(l.get(1));
                    String Chr = hmcntchr.get(chr);
                    if (chr == 43 || chr == 44) {

                    } else {
                        if (chr % 2 == 0) {
                            pos = pos + hmcntsPOS.get(chr);
                        }
                        sb = new StringBuilder();
                        sb.append(Chr).append("\t").append(String.valueOf(pos));
                        bw.write(sb.toString());
                        bw.newLine();
                        if (a % 1000 == 0) {
                            System.out.println("Output " + String.valueOf(a) + " SNPs");
                        }
                        a++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(a) + " SNPs output from " + f.getAbsolutePath());
                System.out.println(chr + " was finished to change the CHROM and Pos from    " + f.getAbsolutePath());
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
     */
    public void extractHapPos(String infileDirS, String outfileDirS) {
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

    /**
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void extractHapPosAllele(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".posAllele.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".posAllele.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnt = 0;
                bw.write("Chr\tPos\tRef\tAlt\n");
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temp = temp.substring(0, 40); //肯定够                 
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4));
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

    /**
     * chr001_5000.vcf --> chr001_5000_IndiHeter.txt Calculate the heterozygote
     * count and propotion by individual taxa
     */
    public void calIndiHeterMiss(String infileDirS, String outfileDirS) {
        //String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/000_sampleVCF";
        //String outfileDirS = "/Users/Aoyue/Documents/test/";
        File[] fs = new File(infileDirS).listFiles();
        //fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            //Start to cal the time beginning
            long startTime = System.nanoTime();
            System.out.println("******************************************************");
            //Start to cal heterozygous
            try {
                String infileS = f.getAbsolutePath();
                //String siteoutfileS = new File(outfileDirS,f.getName().split(".vc")[0] + "_SNPheter.txt").getAbsolutePath();
                String indioutfileS = new File(outfileDirS, f.getName().split(".vc")[0] + "_IndiHeter.txt").getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                //BufferedWriter bw = IOUtils.getTextGzipWriter(siteoutfileS);
                BufferedWriter indibw = IOUtils.getTextWriter(indioutfileS);
                //bw.write("Chr\tPos\tHetNum\tHomNum\tHetPropotion\n");

                String temp;
                String te[] = null;
                /**
                 * *****************定义taxa数组，并添加元素**********************
                 */
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp);
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                }
                /**
                 * *****************定义Genotype
                 * TIntArrayList**********************
                 */
                //TIntArrayList[] genoList = new TIntArrayList[taxa.length];
                //for (int i = 0; i < taxa.length; i++) genoList[i] = new TIntArrayList();
                List[] genoList = new ArrayList[taxa.length];
                for (int i = 0; i < taxa.length; i++) {
                    genoList[i] = new ArrayList();
                }

                while ((temp = br.readLine()) != null) {
                    int genoNum = 0;
                    double homNum = 0;
                    double hetNum = 0;
                    double hetRate = 0;
                    //在一个位点内进行计算

                    if (!temp.startsWith("#")) {
                        //l = PStringUtils.fastSplit(temp, "\t");
                        te = temp.split("\t");
                        if (te[4].length() > 1) {
                            continue;
                        } //only keep biallelic
                        for (int i = 9; i < te.length; i++) {
                            if (te[i].startsWith(".")) {
                                genoList[i - 9].add(0);
                            }
                            if (!te[i].startsWith(".")) {
                                genoNum++; //have the genotype
                                if (te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                    hetNum++; //the number of heterozygous
                                    genoList[i - 9].add(2); // 2 stand for heter
                                }
                                if (te[i].startsWith("0/0") || te[i].startsWith("1/1")) {
                                    homNum++; //the number of heterozygous
                                    genoList[i - 9].add(1); // 1 stand for homo
                                }
                            }
                        }
                        hetRate = hetNum / genoNum;
                        //bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.0f", hetNum)+ "\t"+ String.format("%.0f", homNum) + "\t"+ String.format("%.5f", hetRate) + "\n");
                    }
                }

                /**
                 * ***************** 开始计算个体的杂合度**********************
                 */
                indibw.write("INDV\tTotalSitesWithGeno\tHetSites\tHetProportion\tMissingSites\tMissProportion\n");
                for (int i = 0; i < taxa.length; i++) {
                    double hetSites = Collections.frequency(genoList[i], 2);
                    double totalSitesWithGeno = Collections.frequency(genoList[i], 1) + Collections.frequency(genoList[i], 2); //含有基因型的位点数
                    double hetProportion = hetSites / totalSitesWithGeno;
                    double missSites = Collections.frequency(genoList[i], 0);
                    double missProportion = missSites / (Collections.frequency(genoList[i], 0) + Collections.frequency(genoList[i], 1) + Collections.frequency(genoList[i], 2));
                    indibw.write(taxa[i] + "\t" + String.format("%.0f", totalSitesWithGeno) + "\t" + String.format("%.0f", hetSites) + "\t" + String.format("%.5f", hetProportion) + "\t"
                            + String.format("%.0f", missSites) + "\t" + String.format("%.5f", missProportion) + "\n");
                }
                br.close();
                //bw.flush();
                //bw.close();
                indibw.flush();
                indibw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

            //文件处理完毕，计时
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            //System.out.println("******************************************************" );
            System.out.println(f.getName() + " is finished!!!");
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s" + "    or " + String.format("%.2f", excTime / 60) + " min");
        });
    }

    /**
     * chr001.ABDgenome.filterMiss_subset.vcf.gz
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergesubsetVCF(String infileDirS, String outfileS) {
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        Arrays.sort(fs);
        System.out.println("Chr\tSNP_Num");

        try {
            long startTime = System.nanoTime();
            BufferedReader br = IOUtils.getTextGzipReader(fs[0].getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            br.close();
            for (int i = 0; i < fs.length; i++) {
                br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                int chrint = Integer.parseInt(fs[i].getName().substring(3, 6));
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {

                    } else {
                        cnt++;
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();

                    }
                }
                System.out.println(String.valueOf(chrint) + "\t" + cnt);
            }
            br.close();
            bw.flush();
            bw.close();
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * chr002.ABgenome.filterMiss.vcf --> chr002.ABgenome.filterMiss_subset.vcf
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void subsetVCFRandomParallel(String infileDirS, String outfileDirS, String extractRatio) {

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        //System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        long startTime = System.nanoTime();
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_subset.vcf.gz").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnttotal++;
                        double r = Math.random();
                        double ratio = Double.parseDouble(extractRatio);
                        if (r > ratio) {
                            continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                        }
                        List<String> l = PStringUtils.fastSplit(temp);
                        if (l.get(3).contains(",")) {
                            continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                        }
                        bw.write(temp);
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();

                System.out.println(f.getName() + " with " + cnttotal + " bp has a subset of\t" + cntsubset);
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
        long endTime = System.nanoTime();
        float excTime = (float) (endTime - startTime) / 1000000000;
        System.out.println("Execution time: " + String.format("%.2f", excTime) + "s or  " + String.format("%.2f", excTime / 60) + "min");
        //System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");

    }

    /**
     * chr002.ABgenome.filterMiss.vcf.gz -->
     * chr002.ABgenome.filterMiss_subset.vcf.gz
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void subsetVCFRandomParallel_GZ(String infileDirS, String outfileDirS) {

        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        //System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        long startTime = System.nanoTime();
        fsList.parallelStream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String outfileS = new File(outfileDirS, f.getName().split(".vcf.gz")[0] + "_subset.vcf.gz").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnttotal++;
                        double r = Math.random();
                        if (r > 0.01) {
                            continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                        }
                        List<String> l = PStringUtils.fastSplit(temp);
                        if (l.get(3).contains(",")) {
                            continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                        }
                        bw.write(temp);
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();

                System.out.println(f.getName() + " with " + cnttotal + " bp has a subset of\t" + cntsubset);
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
        long endTime = System.nanoTime();
        float excTime = (float) (endTime - startTime) / 1000000000;
        System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        //System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");

    }

    /**
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void calSNPHetMissMaf(String infileDirS, String outfileDirS) {
        //String infileDirS = "/data4/home/aoyue/vmap2/abd/rawVCF/";
        //String outfileDirS = "/data4/home/aoyue/vmap2/abd/005_vcf/001_calSNPHetMissMaf/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            //开始处理文件，计时
            long startTime = System.nanoTime();
            Calendar cal = Calendar.getInstance();
            int hour = cal.get(Calendar.HOUR_OF_DAY);
            int minute = cal.get(Calendar.MINUTE);
            int second = cal.get(Calendar.SECOND);
            System.out.println("******************************************************");
            System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + "Now starting to " + f.getName() + " heterozygote propotion;");

            try {
                if (f.getName().equals("chr014.vcf")) {

                } else {
                    String infileS = f.getAbsolutePath();
                    String outfileS = new File(outfileDirS, f.getName().split(".vc")[0] + "_SNPheterMiss.txt").getAbsolutePath();

                    BufferedReader br = null;
                    if (infileS.endsWith(".vcf")) {
                        br = IOUtils.getTextReader(infileS);
                    } else if (infileS.endsWith(".vcf.gz")) {
                        br = IOUtils.getTextGzipReader(infileS);
                    }
                    BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                    //bw.write("Chr\tPos\tHetPropotion\n");
                    bw.write("Chr\tPos\tHetNum\tHetPropotion\tMissingNum\tMissProportion\tMaf\n");
                    String temp;
                    String te[] = null;
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
                        //在一个位点内进行计算
                        if (!temp.startsWith("#")) {
                            te = temp.split("\t");
                            if (te[4].length() > 1) {
                                continue;
                            }
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
                                maf = altAF;
                            } else {
                                maf = refAF;
                            }
                            //bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.5f", hetRate) + "\n");
                            bw.write(te[0] + "\t" + te[1] + "\t" + String.format("%.0f", hetNum) + "\t" + String.format("%.5f", hetRate)
                                    + "\t" + String.format("%.0f", missNum) + "\t" + String.format("%.5f", missRate) + "\t" + String.format("%.5f", maf) + "\n");
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                }

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

            //文件处理完毕，计时
            hour = cal.get(Calendar.HOUR_OF_DAY);
            minute = cal.get(Calendar.MINUTE);
            second = cal.get(Calendar.SECOND);
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            //System.out.println("******************************************************" );
            //System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + f.getName() + " is finished!!!");
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        });

    }

    /**
     *
     * @param infileS
     * @param outfileS
     */
    private void subsetVCFdataRandom_singleStream(String infileS, String outfileS) {
        //String infileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001.vcf";
        //String outfileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001_subset.vcf";

        //String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/002_merge/subsetchr1_15.vcf.gz";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/002_merge/subset10ksnp.vcf.gz";
        try {
            //BufferedReader br = IOUtils.getTextReader(infileS);
            //BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);

            String temp = null;
            int cnt = 0;
            System.out.println(new SimpleDateFormat().format(new Date()) + "    program execution.\n");
            long startTime = System.nanoTime();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                } else {
                    cnt++;
                    double r = Math.random();
                    //if (r > 0.001) {
                    if (r > 0.067) {
                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).contains(",")) {
                        continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                    }
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");

            System.out.println("Chr 1 snp number is " + cnt + ".");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * chr005.Dgenome.vcf.gz --> chr005.Dgenome.bi.vcf.gz
     *
     * @param infileDirS
     */
    public void filterAllele(String infileDirS, String outfileDirS) {
        new File(outfileDirS).mkdirs();
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println("Chr\tSNPNum\tBiallelicNum\tIndelNum\tInsertionNum\tDelectionNum\t");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".bi.vcf.gz")).getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnt = 0;
                int snpNum = 0;
                int biallelicNum = 0;
                int indelNum = 0;
                int insertionNum = 0;
                int delectionNum = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnt++;
                        String alt = PStringUtils.fastSplit(temp).get(4);
                        if (alt.contains("D")) {
                            delectionNum++;
                        }
                        if (alt.contains("I")) {
                            insertionNum++;
                        }

                        if (!(alt.contains(",")) && !(alt.equals("D")) && !(alt.equals("I"))) {
                            biallelicNum++;
                            bw.write(temp);
                            bw.newLine();
                        }

                    }

                }
                indelNum = delectionNum + insertionNum;
                snpNum = cnt - indelNum;
                br.close();
                bw.flush();
                bw.close();
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(snpNum) + "\t" + String.valueOf(biallelicNum) + "\t" + String.valueOf(indelNum) + "\t" + String.valueOf(insertionNum) + "\t" + String.valueOf(delectionNum));
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 对已生成的14条染色体进行计数.注意加log文件，结果在log文件中显示
     *
     * @param infileDirS
     */
    public void countSitesinFastCallformat(String infileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);

        System.out.println("Chr\tSNPNum\tBiallelicNum\tIndelNum\tInsertionNum\tDelectionNum\t");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = null;
                int cnt = 0;
                int snpNum = 0;
                int biallelicNum = 0;
                int indelNum = 0;
                int insertionNum = 0;
                int delectionNum = 0;
                while ((temp = br.readLine()) != null) { //是否含有D I
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    cnt++;
                    String alt = PStringUtils.fastSplit(temp).get(4);
                    if (alt.contains("D")) {
                        delectionNum++;
                    }
                    if (alt.contains("I")) {
                        insertionNum++;
                    }

                    if (!(alt.contains(",")) && !(alt.equals("D")) && !(alt.equals("I"))) {
                        biallelicNum++;
                    }
                }

                indelNum = delectionNum + insertionNum;
                snpNum = cnt - indelNum;
                br.close();
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(snpNum) + "\t" + String.valueOf(biallelicNum) + "\t" + String.valueOf(indelNum) + "\t" + String.valueOf(insertionNum) + "\t" + String.valueOf(delectionNum));
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    /**
     * 将计算出的snp位点数进行合并，成1D 2D 3D 4D 5D 6D 7D形式；
     */
    public void mergeChr1and2_Dgenome(String infileS, String outfileS) {
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites_mergeChr1and2.txt";

        //infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK.txt";
        //outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK_mergeChr1and2.txt";
        String[] chr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
        HashMap<Integer, String> hmcntchr = new HashMap<>();
        HashMap<Integer, Integer> hmcntSNPNum = new HashMap<>();
        for (int i = 0; i < chr.length; i++) {
            hmcntchr.put(cnts[i], chr[i]);
        }

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine(); //read header
            bw.write(temp);
            bw.newLine(); //writer header
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                int site1 = Integer.parseInt(PStringUtils.fastSplit(temp).get(1));
                if ((temp = br.readLine()) != null) {
                    int site2 = Integer.parseInt(PStringUtils.fastSplit(temp).get(1));
                    int site = site1 + site2;
                    cnt++;
                    hmcntSNPNum.put(cnt, site);
                    bw.write(hmcntchr.get(cnt) + "\t" + hmcntSNPNum.get(cnt));
                    bw.newLine();
                } else {

                }
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
     * 将计算出的snp位点数进行合并，成1A 1B 1D形式；
     */
    public void mergeChr1and2(String infileS, String outfileS) {
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites_mergeChr1and2.txt";

        infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK.txt";
        outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK_mergeChr1and2.txt";
        String[] chr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
        HashMap<Integer, String> hmcntchr = new HashMap<>();
        HashMap<Integer, Integer> hmcntSNPNum = new HashMap<>();
        for (int i = 0; i < chr.length; i++) {
            hmcntchr.put(cnts[i], chr[i]);
        }

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine(); //read header
            bw.write(temp);
            bw.newLine(); //writer header
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                int site1 = Integer.parseInt(PStringUtils.fastSplit(temp).get(1));
                if ((temp = br.readLine()) != null) {
                    int site2 = Integer.parseInt(PStringUtils.fastSplit(temp).get(1));
                    int site = site1 + site2;
                    cnt++;
                    hmcntSNPNum.put(cnt, site);
                    bw.write(hmcntchr.get(cnt) + "\t" + hmcntSNPNum.get(cnt));
                    bw.newLine();
                } else {

                }
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
     * Count the snp sites via parallel stream and print it into the
     * inputstream.
     *
     * @param infileDirS Usage:java -jar PlantGenetics.jar > countSites.txt &
     */
    public void countSites_parallelStream(String infileDirS) {
        //infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/002_subsetVCF/001_subsetVCF/";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        //System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        System.out.println("Chr\tSNPNum");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                String infileS = f.getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    cnt++;
                }
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cnt));
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        //System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");
    }

    /**
     * Count the snp sites via single stream and print it into the inputstream.
     *
     * @param infileDirS
     */
    public void countSites_singleStream(String infileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    cnt++;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(String.valueOf(cnt) + "\t" + fs[i].getName());
            sum += cnt;
        }
        System.out.println(String.valueOf(sum));
    }
}
