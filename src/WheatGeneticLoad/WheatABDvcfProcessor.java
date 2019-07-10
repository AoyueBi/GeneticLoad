/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
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
public class WheatABDvcfProcessor {

    public WheatABDvcfProcessor() {
        //this.subsetVCFdataRandom();
        //new Treetest();
        //this.countSites("/data4/home/aoyue/vmap2/abd/rawVCF/");
        //this.calSNPHetMissMaf();
        this.mergeChr1and2();

    }

    public void mergeChr1and2() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites_mergeChr1and2.txt";
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

    public void calSNPHetMissMaf() {
        String infileDirS = "/data4/home/aoyue/vmap2/abd/rawVCF/";
        String outfileDirS = "/data4/home/aoyue/vmap2/abd/005_vcf/001_calSNPHetMissMaf/";
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
     * 对已生成的12条染色体进行计数. CMD: java -jar PlantGenetics.jar > countSites.txt & [2]
     * 350108
     *
     * @param infileDirS
     */
    public void countSites(String infileDirS) {
        //infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        System.out.println("Chr\tSNP Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().split("chr")[1].split(".vcf")[0]; //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    cnt++;
                }
                br.close();
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cnt));
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");
    }

    private void subsetVCFdataRandom() {
        String infileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001.vcf";
        String outfileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001_subset.vcf";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
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
                    if (r > 0.001) {
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
}
