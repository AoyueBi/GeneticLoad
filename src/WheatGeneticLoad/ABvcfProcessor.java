/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.CountSites;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;

/**
 *
 * @author Aoyue
 */
public class ABvcfProcessor {

    public ABvcfProcessor() {
//        new CountSites().filterSNPtoBi("/data4/home/aoyue/vmap2/analysis/001_rawvcf/ab/", "/data4/home/aoyue/vmap2/analysis/002_bivcf/ab/");
//        this.subsetVCFRandomParallel();
//        this.mergesubsetVCF();
//        new Treetest().labels();
//        new Treetest().TREE_COLORS();
//        new Treetest().colRangesbyDico();
        this.addGrouptoMDS();
        //this.checkIBSdistance();

        this.getAllValuefromMatrix();
        //this.removeRepeat();
        //this.getsubsetMatrix();
        this.getValuefrom125IBSdistance();

    }

    /**
     * 从125matrix得到两两之间关系矩阵的值，列成一列，为了画密度分布图
     */
    public void getValuefrom125IBSdistance() {
        String distancefileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/007_getsubset125IBSdistance.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/008_IBSdistance_S125.txt";
        int taxanum = 125;
        try {
            BufferedReader br = IOUtils.getTextReader(distancefileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();

            String[][] matrix = new String[125][125];
            int r = -1;
            //建立matrix
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                r++; //r代表了矩阵行的index； j 代表了第 
                for (int j = 0; j < taxanum; j++) {
                    matrix[r][j] = l.get(j + 1);
                }
            }
            br.close();

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("IBSdistance");
            bw.newLine();
            for (int i = 0; i < taxanum; i++) {
                for (int j = 0; j < i + 1; j++) {
                    bw.write(matrix[i][j]);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 从205个样品的matrix中，获取不是重复的额165个样品两两比较的matrix，然后将该matrix种的值提出出来方成一类。
     */
    public void getsubsetMatrix() {
        String repeatlistS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_labelsChange/reaptItemlist.txt";
        String distancefileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/006_matrix_S205.modify.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/007_getsubset125IBSdistance.txt";
        List<String> rep = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(repeatlistS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                rep.add(temp);
            }
            String[] reps = rep.toArray(new String[rep.size()]);
            Arrays.sort(reps);

            //将原始matrix读进表格
            RowTable<String> t = new RowTable<>(distancefileS);
            boolean[] ifout = new boolean[t.getRowNumber()];
            List<Integer> indexColum = new ArrayList<>();

            //开始删除列
            for (int i = 0; i < reps.length; i++) {
                t.removeColumn(reps[i]);
            }

            //开始删除行
            for (int i = 0; i < t.getRowNumber(); i++) {
                String rowname = t.getCell(i, 0);
                int f = Arrays.binarySearch(reps, rowname);
                if (f < 0) {
                    ifout[i] = true;
                }
            }

            t.writeTextTable(outfileS, IOFileFormat.Text, ifout);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 除去40对重复的IBS distance值，但是这种想法是错误的。呵呵呵额
     */
    public void removeRepeat() {
        String allIBSfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/004_IBSdistance_S205.txt";
        String repeatfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/001_checkIBSdistance.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/005_removeRepeatIBSdistance.txt";
        try {
            RowTable<String> t = new RowTable<>(allIBSfileS);
            List<String> db = t.getColumn(5);
            Collections.sort(db);

            BufferedReader br = IOUtils.getTextReader(allIBSfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                int index = Collections.binarySearch(db, temp);
                if (index < 0) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 获取matrix中所有两两比较的遗传距离值，放入一列中，为了做密度分布图用。
     */
    public void getAllValuefromMatrix() {
        int taxanum = 0;
        String distancefileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/001_newick/chr1_4.ABgenome.filterMiss_subset1.8m.matrix.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/004_IBSdistance_S205.txt";

        try {
            BufferedReader br = IOUtils.getTextReader(distancefileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            //确定taxa数量
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                taxanum = Integer.parseInt(temp);
                break;
            }

            String[][] matrix = new String[taxanum][taxanum];
            int r = -1;
            //建立matrix
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                r++;
                for (int j = 0; j < taxanum; j++) {
                    matrix[r][j] = l.get(j + 1);
                }
            }
            br.close();

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("IBSdistance");
            bw.newLine();
            for (int i = 0; i < taxanum; i++) {
                for (int j = 0; j < i + 1; j++) {
                    bw.write(matrix[i][j]);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 获取重复的40个样品，两两比较的IBS distance值
     *
     * ##IBS_Distance_Matrix.AverageTotalSites=1406356.8265901483
     * ##IBS_Distance_Matrix.NumAlleles=3 ##IBS_Distance_Matrix.TrueIBS=false
     * ##Matrix_Type=IBS_Distance_Matrix 
     * ##205 
     * ##B0001
     */
    public void checkIBSdistance() {
        int taxanum = 0;
        String distancefileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/001_newick/chr1_4.ABgenome.filterMiss_subset1.8m.matrix.txt";
        String reaptItemS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/reaptItem.txt";
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_AB_S205_germplasmInfo.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/ab/002_checkIBSdistance/001_checkIBSdistance.txt";
        HashMap<String, String> hmDatabaseIDContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDPartContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDCty = new HashMap<>();
        HashMap<String, Integer> hmContinentIndex = new HashMap<>();
        HashMap<String, String> hmEmmertype = new HashMap<>();
        HashMap<String, String> hmtaxa = new HashMap<>();

        try {
            //建立数组
            BufferedReader br = IOUtils.getTextReader(distancefileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                taxanum = Integer.parseInt(temp); //matrix文件显示taxa数目
                break;
            }
            String[][] matrix = new String[taxanum][taxanum];
            int r = -1;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                r++; //r从0开始，j也从0开始，matrix的二维矩阵就是00
                for (int j = 0; j < taxanum; j++) {
                    matrix[r][j] = l.get(j + 1);
                }
            }
            br.close();

            //建立种质信息数组
            RowTable<String> t = new RowTable<>(dbfileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                String databaseID = t.getCell(i, 0);
                String continent = t.getCell(i, 12);
                String partContinent = t.getCell(i, 14);
                String country = t.getCell(i, 10);
                String type = t.getCell(i, 7);
                String taxa = t.getCell(i, 4);
                hmDatabaseIDContinent.put(databaseID, continent);
                hmDatabaseIDPartContinent.put(databaseID, partContinent);
                hmDatabaseIDCty.put(databaseID, country);
                hmEmmertype.put(databaseID, type);
                hmtaxa.put(databaseID, taxa);
            }

            br = IOUtils.getTextReader(reaptItemS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tDatabase1\tDatebase2\tType\tCty\tIBSdistance");
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                String dbid1 = PStringUtils.fastSplit(temp).get(1).substring(0, 5);
                String dbid2 = PStringUtils.fastSplit(temp).get(1).substring(6);
                bw.write(hmtaxa.get(dbid1) + "\t" + dbid1 + "\t" + dbid2 + "\t" + hmEmmertype.get(dbid1) + "\t" + hmDatabaseIDCty.get(dbid1) + "\t" + matrix[Integer.parseInt(dbid1.substring(1)) - 1][Integer.parseInt(dbid2.substring(1)) - 1]);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    // 对MDS方法生成的PC结果添加分组信息，在第二列中加入 
    //String[] group = {"Oceania","Africa","North America","South America","Europe","Central Asia","South Asia","Western Asia","East Asia","NA"};
    //String[] col = {"#F1E1FF","#F4D03F","#F1948A","#5DADE2","#ABEBC6","#239B56","#CD6155","#FF6347","#7B241C","#EBEDEF"};
    //String[] col = {"#5DADE2","#7B241C","#F1E1FF","#F4D03F","#FF9900","#006600","#389038","#82C782","#CCFFCC","#EBEDEF"};
    //分组1为：indexColum 
    //分组2为：大洋洲 非洲 北美洲 南美洲 欧洲 亚洲 
    //"Oceania","Africa","North America","South America","Europe","Asia" 
    //颜色为："#F1E1FF","#F4D03F","#5DADE2","#7B241C","#FF9900","#82C782"
    // 分组3为：大洋洲 非洲 北美洲 南美洲 欧洲部洲 亚洲部洲
    //分组4为国家
    //0 <- "Africa" 
    // 1 <- "Asia" 
    // 2 <- "Europe" 
    // 3 <- "North America" 
    // 4 <- * "Oceania"
    // 5 <- "South America"
    //6<- "NA"
    public void addGrouptoMDS() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/005_pca/ab/001_mdsMethod/002_MDS_PCs_Matrix_subset63ksnp_forR.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/005_pca/ab/001_mdsMethod/003_MDS_PCs_Matrix_subset63ksnp_forR_addGroup.txt";
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_AB_S205_germplasmInfo.txt";
        HashMap<String, String> hmDatabaseIDContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDPartContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDCty = new HashMap<>();
        HashMap<String, Integer> hmContinentIndex = new HashMap<>();
        HashMap<String, String> hmEmmertype = new HashMap<>();

        RowTable<String> t = new RowTable<>(dbfileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String databaseID = t.getCell(i, 0);
            String continent = t.getCell(i, 12);
            String partContinent = t.getCell(i, 14);
            String country = t.getCell(i, 10);
            String type = t.getCell(i, 7);
            hmDatabaseIDContinent.put(databaseID, continent);
            hmDatabaseIDPartContinent.put(databaseID, partContinent);
            hmDatabaseIDCty.put(databaseID, country);
            hmEmmertype.put(databaseID, type);
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
            bw.write("GroupIndex\tCountry\tContinent\tPart-Cpntinent\tType");
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
                    bw.write("1\tCHN\tAsia\tEast Asia\tdicoccum");
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + l.get(i));
                    }
                    bw.newLine();

                } else {
                    if (cty.equals("NA")) {
                        bw.write(l.get(0) + "\t");
                        bw.write("6\tNA\tNA\tNA\t" + hmEmmertype.get(l.get(0)));
                        for (int i = 1; i < l.size(); i++) {
                            bw.write("\t" + l.get(i));
                        }
                        bw.newLine();

                    } else {
                        bw.write(id + "\t");
                        bw.write(String.valueOf(hmContinentIndex.get(hmDatabaseIDContinent.get(id))) + "\t" + hmDatabaseIDCty.get(id) + "\t" + hmDatabaseIDContinent.get(id)
                                + "\t" + hmDatabaseIDPartContinent.get(id) + "\t" + hmEmmertype.get(l.get(0)));
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
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/002_subsetVCF/001_subsetVCF";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/002_subsetVCF/002_merge/chr1_4.ABgenome_subset.vcf.gz";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        Arrays.sort(fs);
        System.out.println("Chr\tSNP_Num");

        try {
            BufferedReader br = IOUtils.getTextGzipReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/002_subsetVCF/001_subsetVCF/chr001.ABgenome_subset.vcf.gz");
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

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        /**
         * Chr	SNP_Num 1	20890 2	7296 3	20311 4	14748
         */
    }

    /**
     * 对已生成的vcf数据进行随机抽取，多线程n条染色体同时进行；
     */
    public void subsetVCFRandomParallel() {
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/003_filterMiss/ab/";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/004_subsetvcf/ab/";

        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        long startTime = System.nanoTime();
        fsList.parallelStream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String outfileS = new File(outfileDirS, f.getName().substring(0, 15) + "_subset.vcf").getAbsolutePath();
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

}
