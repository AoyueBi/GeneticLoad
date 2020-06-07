/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.CountSites;
import AoUtils.Script;
import pgl.infra.table.RowTable;
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
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class ABDvcfProcessor {

    public ABDvcfProcessor() {
//        this.subsetVCFdataRandom();
//        new Treetest();
//        this.countSites("/data4/home/aoyue/vmap2/abd/rawVCF/");
//        this.calSNPHetMissMaf();
        //this.mergeChr1and2();
//        this.subsetHetMissMafRandom(); //已和下一步合并起来
//        this.mergesubsetHetMissMafRandomData();
//
//        this.subsetVCFRandomParallel();
        //this.mergesubsetVCF();
//        
        //this.addGrouptoMDS();
//        this.VCFfromGATKTest();
//        new Script().bgzip_noscript("/data4/home/aoyue/vmap2/abd/rawVCF/", "/data4/home/aoyue/vmap2/abd/rawVCF/");
//        new CountSites().countSitesinFastCallformat("/data4/home/aoyue/vmap2/analysis/001_rawvcf/abd");
//        new CountSites().filterSNPtoBi("/data4/home/aoyue/vmap2/analysis/001_rawvcf/abd/", "/data4/home/aoyue/vmap2/analysis/002_bivcf/abd/");
        //new CountSites().subsetVCFRandomParallel_GZ("/data4/home/aoyue/vmap2/analysis/003_filterMiss/abd/", "/data4/home/aoyue/vmap2/analysis/004_subsetvcf/abd/");
        
        
    }
    
    /**
     * 对CATK生成的VCF文件进行测试，看看结果是否和FastCall一样
     */
    public void VCFfromGATKTest(){
        this.calSNPHetMissMaffromGATK();
        this.countSitesfromGATK("/data3/wgs/vcf/GATK/001_rawvcf/SNP_S319");
        
    }
    
    /**
     * 对已生成的42条染色体进行计数. CMD: java -jar PlantGenetics.jar > countSites_fromGATK.txt & 
     * [1] 59847
     *
     * @param infileDirS
     */
    public void countSitesfromGATK(String infileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        System.out.println("Chr\tSNP Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().split(".snp.vcf")[0]; //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
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
    
    public void calSNPHetMissMaffromGATK() {
        String infileDirS = "/data3/wgs/vcf/GATK/001_rawvcf/SNP_S319/";
        String outfileDirS = "/data4/home/aoyue/vmap2/abd/GATK/001_calSNPHetMissMaf/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
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
                    String infileS = f.getAbsolutePath();
                    String outfileS = new File(outfileDirS, f.getName().split(".snp.vc")[0] + "_SNPheterMiss.txt.gz").getAbsolutePath();

                    BufferedReader br = null;
                    if (infileS.endsWith(".vcf")) {
                        br = IOUtils.getTextReader(infileS);
                    } else if (infileS.endsWith(".vcf.gz")) {
                        br = IOUtils.getTextGzipReader(infileS);
                    }
                    BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
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
     * 对MDS方法生成的PC结果添加分组信息，在第二列中加入 String[] group = {"Oceania","Africa","North
     * America","South America","Europe","Central Asia","South Asia","Western
     * Asia","East Asia","NA"}; //String[] col =
     * {"#F1E1FF","#F4D03F","#F1948A","#5DADE2","#ABEBC6","#239B56","#CD6155","#FF6347","#7B241C","#EBEDEF"};
     * String[] col =
     * {"#5DADE2","#7B241C","#F1E1FF","#F4D03F","#FF9900","#006600","#389038","#82C782","#CCFFCC","#EBEDEF"};
     * 分组1为：index 分组2为：大洋洲 非洲 北美洲 南美洲 欧洲 亚洲 "Oceania","Africa","North
     * America","South America","Europe","Asia"
     * 颜色为："#F1E1FF","#F4D03F",    "#5DADE2","#7B241C","#FF9900","#82C782" 分组3为：大洋洲
     * 非洲 北美洲 南美洲 欧洲部洲 亚洲部洲 分组4为国家
     *
     * 0 <- "Africa" 1 <- "Asia" 2 <- "Europe" 3 <- "North America" 4 <-
     * "Oceania" 5 <- "South America" 6<- "NA"
     */
    public void addGrouptoMDS() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/003_mdsMethod/002_MDS_PCs_Matrix_subset10ksnp_forR.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/003_mdsMethod/003_MDS_PCs_Matrix_subset10ksnp_forR_addGroup.txt";
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/002_merge/002_All373wheat_ABD_CountryBreedingStatus_Eng.txt";
        HashMap<String, String> hmDatabaseIDContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDPartContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDCty = new HashMap<>();
        HashMap<String, Integer> hmContinentIndex = new HashMap<>();

        RowTable<String> t = new RowTable<>(dbfileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String databaseID = t.getCell(i, 2);
            String continent = t.getCell(i, 6);
            String partContinent = t.getCell(i, 8);
            String country = t.getCell(i, 4);
            hmDatabaseIDContinent.put(databaseID, continent);
            hmDatabaseIDPartContinent.put(databaseID, partContinent);
            hmDatabaseIDCty.put(databaseID, country);
        }

        String[] continents = {"Africa", "Asia", "Europe", "North America", "Oceania", "South America"};
        for (int i = 0; i < continents.length; i++) {
            hmContinentIndex.put(continents[i], i);
        }

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine(); //header
            List<String> l = PStringUtils.fastSplit(temp);
            bw.write(l.get(0) + "\t");
            bw.write("GroupIndex\tCountry\tContinent\tPart-Cpntinent");
            for (int i = 1; i < l.size(); i++) {
                bw.write("\t" + l.get(i));
            }
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String id = l.get(0);
                String cty = hmDatabaseIDCty.get(id);
                if (id.equals("CS")) {
                    bw.write(l.get(0) + "\t");
                    bw.write("1\tCHN\tAsia\tEast Asia");
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + l.get(i));
                    }
                    bw.newLine();

                } else {
                    if (cty.equals("NA")) {
                        bw.write(l.get(0) + "\t");
                        bw.write("6\tNA\tNA\tNA");
                        for (int i = 1; i < l.size(); i++) {
                            bw.write("\t" + l.get(i));
                        }
                        bw.newLine();

                    } else {
                        bw.write(id + "\t");
                        bw.write(String.valueOf(hmContinentIndex.get(hmDatabaseIDContinent.get(id))) + "\t" + hmDatabaseIDCty.get(id) + "\t" + hmDatabaseIDContinent.get(id)
                                + "\t" + hmDatabaseIDPartContinent.get(id));
                        for (int i = 1; i < l.size(); i++) {
                            bw.write("\t" + l.get(i));
                        }
                        bw.newLine();

                    }

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
     * 对抽取的样本进行合并，并统计每个文件抽取了多少条SNP
     */
    public void mergesubsetVCF() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/001_subsetVCF/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/002_merge/subsetchr1_15.vcf.gz";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        Arrays.sort(fs);
        System.out.println("Chr\tSNP_Num");

        try {
            BufferedReader br = IOUtils.getTextGzipReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/001_subsetVCF/chr001_subset.vcf.gz");
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
                int chrint = Integer.parseInt(fs[i].getName().split("chr")[1].split("_subset")[0]);
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

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        /**
         * Chr	SNP_Num 1	12291 2	4260 3	14307 4	8947 5	11613 6	1231 7	11787 8
         * 9473 9	15338 10	12768 11	11306 12	5451 13	9547 14	9149 15	13244 Total
         * 150712
         */

    }

    /**
     * 对已生成的vcf数据进行随机抽取，多线程14条染色体同时进行；
     */
    public void subsetVCFRandomParallel() {
        //String infileDirS = "/data4/home/aoyue/vmap2/abd/rawVCF/";
        //String outfileDirS = "/data4/home/aoyue/vmap2/abd/005_vcf/004_pca/001_subsetVCF/";
        
        String infileDirS = "";
        String outfileDirS = "";
        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        long startTime = System.nanoTime();
        fsList.parallelStream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_subset.vcf").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        
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
                        cnt++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();

                System.out.println(f.getName() + " is being subset about\t" + cnt);
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
        long endTime = System.nanoTime();
        float excTime = (float) (endTime - startTime) / 1000000000;
        System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");

    }

    /**
     * 将生成的1-6染色体统计数据合并起来，成为一个文件。 1.建立表头，写入头文件； 2.建立数组，将文件名字按顺序读入；
     * 3.进行for循环，合并文件。
     */
    public void mergesubsetHetMissMafRandomData() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_calSNPHetMissMaf/002_subsetCalSNPHetMissMaf/";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_calSNPHetMissMaf/003_merge/subsetCalSNPHetMissMaf.txt";
        
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/003_calSNPHetMissMaf/002_subsetCalSNPHetMissMaf/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/003_calSNPHetMissMaf/003_merge/subsetCalSNPHetMissMaf.txt";

        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        Arrays.sort(fs);

        try {
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_calSNPHetMissMaf/002_subsetCalSNPHetMissMaf/chr001_SNPheterMiss_subset.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            br.close();
            //按文件顺序合并
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                br = IOUtils.getTextReader(infileS);
                String temp = br.readLine(); //表头
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
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
     * 从全部数据中提取Chr HetPropotion MissProportion Maf 4列，合计1M数据量
     */
    public void subsetHetMissMafRandom() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_calSNPHetMissMaf";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_subsetCalSNPHetMissMaf";
        
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/003_calSNPHetMissMaf/001_allpos/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/003_calSNPHetMissMaf/002_subsetCalSNPHetMissMaf/";

        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        fsList.parallelStream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_subset.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = br.readLine();
                bw.write("Chr\tHetPropotion\tMissProportion\tMaf\n");
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    double r = Math.random();
                    if (r > 0.001) {
                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    bw.write(l.get(0) + "\t" + l.get(3) + "\t" + l.get(5) + "\t" + l.get(6));
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");
        this.mergesubsetHetMissMafRandomData();

    }

    /**
     * 将计算出的snp位点数进行合并，成1A 1B 1D形式；
     */
    public void mergeChr1and2() {
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites_mergeChr1and2.txt";
        
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK_mergeChr1and2.txt";
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
}
