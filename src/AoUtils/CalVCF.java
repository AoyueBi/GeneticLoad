/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.table.TableInterface;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class CalVCF {

    public CalVCF() {
//        this.calSiteMeanDepth();
//        this.reheader();
//        this.mkPhylipFormat();
//        this.mkdirs();
//        this.extractVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF/chr032_exon_vmap2.1.vcf.gz","/Users/Aoyue/Documents/test.vcf","/Users/Aoyue/Documents/test.txt");

//        this.getIBSdistance();

    }

    /**
     * chr001.ABDgenome.filterMiss_subset.vcf.gz
     *
     * @param infileDirS
     * @param outfileS
     */
    public static void mergeVCF(String infileDirS, String outfileS) {
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            long startTime = System.nanoTime();
            BufferedReader br = AoFile.readFile(fs[0].getAbsolutePath());
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            br.close();
            int total = 0;
            for (int i = 0; i < fs.length; i++) {
                br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
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
                total = total + cnt;
                System.out.println(cnt + "\tsnps in " + fs[i].getName());
            }
            System.out.println(total + "\tsnps totally, mergevcf pipeline is completed at\t" + outfileS);
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
     *
     */
    public void getIBSdistance(){

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/Dsubgenome_diploid.vcf.gz";
        GenotypeGrid gt = new GenotypeGrid(infileS,GenoIOFormat.VCF_GZ);

        int taxaNum = gt.getTaxaNumber();
        String[] taxa = gt.getTaxaNames();
        System.out.println(taxaNum + "\t" + taxa.length);
        float[][] ibsMatrix =   gt.getIBSDistanceMatrix();

//        GenotypeGrid[] gts = new GenotypeGrid[fN];
//        int totalSiteCount = 0;
//        for (int i = 0; i < gts.length; i++) { //对genotypeTable 进行初始化
//            gts[i] = new GenotypeGrid(fList.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
//            totalSiteCount+=gts[i].getSiteNumber(); //获取 fN 个文件的总位点数
//        }
//        GenotypeOperation.mergeGenotypesBySite(gts[0], gts[1]); // 合并两个VCF文件
//        GenotypeGrid gt = gts[0];
//        TDoubleArrayList missingSite = new TDoubleArrayList(); // calculation 1
//        TDoubleArrayList hetSite = new TDoubleArrayList(); // calculation 2
//        TDoubleArrayList maf = new TDoubleArrayList(); // calculation 3
//        TDoubleArrayList missingTaxon = new TDoubleArrayList(); // calculation 4
//        TDoubleArrayList hetTaxon = new TDoubleArrayList(); // calculation 5
//
//        int step = gt.getSiteNumber()/size;
//        for (int i = 0; i < gt.getSiteNumber(); i+=step) {
//            missingSite.add(((double) gt.getMissingNumberBySite(i)/gt.getTaxaNumber()));
//            hetSite.add(gt.getHeterozygousProportionBySite(i));
//            maf.add(gt.getMinorAlleleFrequency(i));
//        }
//        for (int i = 0; i < gt.getTaxaNumber(); i++) {
//            missingTaxon.add((double)gt.getMissingNumberByTaxon(i)/gt.getSiteNumber());
//            hetTaxon.add(gt.getHeterozygousProportionByTaxon(i));
//        }
//
//        String siteQCfileS = new File(outDirS,genomeType + "_site_QC.txt.gz").getAbsolutePath();
//        String taxaQCFileS = new File (outDirS, genomeType+"_taxa_QC.txt.gz").getAbsolutePath();

    }

    /**
     *
     * @param infileS
     * @param taxaArray
     * @return
     */
    public static List<String> extractGenotable(String infileS, List<String> taxaArray){
        List<String> out = new ArrayList<>();
        List<Integer> indexTaxa = new ArrayList<>();
        Collections.sort(taxaArray);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    String element = "CHROM" + "\t" + l.get(1);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 2; i < l.size(); i++) { //需要修改需要修改需要修改需要修改需要修改需要修改需要修改
                        String taxon = l.get(i);
                        int index1 = Collections.binarySearch(taxaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                            sb.append("\t").append(l.get(i));
                        }
                    }
                    String line = element + sb.toString();
                    out.add(line);
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("CHROM")) {
                    l = PStringUtils.fastSplit(temp);

                    String chr = l.get(0);
                    String pos = l.get(1);
                    String element = chr + "\t" + pos;
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        String geno = l.get(indexTaxa.get(i));
                        sb.append("\t").append(geno);
                    }
                    String line = element + sb.toString();
                    out.add(line);
                }
            }
            br.close();

            System.out.println(infileS + " is completed with line number (with header)" + out.size() + "\tActual taxa size: " + indexTaxa.size() + "\tTotal sites : " );
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }

    public static void extractGenotable(String infileS, List<String> taxaArray, String outfileS){

        List<Integer> indexTaxa = new ArrayList<>();
        Collections.sort(taxaArray);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write("CHROM" + "\t" + l.get(1));
                    for (int i = 2; i < l.size(); i++) { //需要修改需要修改需要修改需要修改需要修改需要修改需要修改
                        String taxon = l.get(i);
                        int index1 = Collections.binarySearch(taxaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                            bw.write("\t" + l.get(i));
                        }
                    }
                    bw.newLine(); //写完之后记得换行
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("CHROM")) {
                    l = PStringUtils.fastSplit(temp);

                    String chr = l.get(0);
                    String pos = l.get(1);
                    bw.write(chr + "\t" + pos);

                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        String geno = l.get(indexTaxa.get(i));
                        bw.write("\t" + geno);
                    }
                    bw.newLine();
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexTaxa.size() + "\tTotal sites : " + cnt);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     *
     * subset genotype table from total genotype table
     *
     * @param infileS VCF file
     * @param taxaArray taxalist without header
     * @param outfileS TXT file
     */
    public static void extractGenotable(String infileS, String[] taxaArray, String outfileS){

        List<Integer> indexTaxa = new ArrayList<>();
        Arrays.sort(taxaArray);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write("CHROM" + "\t" + l.get(1));
                    for (int i = 2; i < l.size(); i++) { //需要修改需要修改需要修改需要修改需要修改需要修改需要修改
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(taxaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                            bw.write("\t" + l.get(i));
                        }
                    }
                    bw.newLine(); //写完之后记得换行
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("CHROM")) {
                    l = PStringUtils.fastSplit(temp);

                    String chr = l.get(0);
                    String pos = l.get(1);
                    bw.write(chr + "\t" + pos);

                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        String geno = l.get(indexTaxa.get(i));
                        bw.write("\t" + geno);
                    }
                    bw.newLine();
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexTaxa.size() + "\tTotal sites : " + cnt);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * convert specific taxa in taxaListFileS(without header) into genotype table,
     * aims to calculate the heterozygous genotype proportion along chromosome later or for other pipeline.
     * Here, I define 0/0 -> 0, 0/1 -> 1, 1/1 -> 2, ./. -> NA
     * @param infileS VCF file
     * @param taxaArray taxalist without header
     * @param outfileS TXT file
     */
    public static void extractVCFtable(String infileS, String[] taxaArray, String outfileS){

        List<Integer> indexTaxa = new ArrayList<>();
        Arrays.sort(taxaArray);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中

                }
                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write("CHROM" + "\t" + l.get(1));
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(taxaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                            bw.write("\t" + l.get(i));
                        }
                    }
                    bw.newLine(); //写完之后记得换行
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    String altList = l.get(4);
                    int nAlt = PStringUtils.fastSplit(altList, ",").size();
                    if (nAlt > 1) continue; //filter alt num with 2 or more
                    if (altList.equals("D") || altList.equals("I")) continue; //filter D I

                    String chr = l.get(0);
                    String pos = l.get(1);
                    bw.write(chr + "\t" + pos);

                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        String geno = l.get(indexTaxa.get(i));
                        geno = PStringUtils.fastSplit(geno,":").get(0);
//                        System.out.println(geno);
                        if(geno.equals("0/0")){
                            bw.write( "\t0");
                        }
                        if(geno.equals("0/1") || geno.equals("1/0")){
                            bw.write(  "\t1");
                        }
                        if(geno.equals("1/1")){
                            bw.write( "\t2");
                        }
                        if(geno.equals("./.")){
                            bw.write("\t9");
                        }
                    }
                    bw.newLine();
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexTaxa.size() + "\tTotal sites : " + cnt);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * convert specific taxa in taxaListFileS(without header) into genotype table,
     * aims to calculate the heterozygous genotype proportion along chromosome later or for other pipeline.
     * Here, I define 0/0 -> 0, 0/1 -> 1, 1/1 -> 2, ./. -> NA
     * @param infileS VCF file
     * @param taxaListFileS taxalist without header
     * @param outfileS TXT file
     */
    public static void extractVCFtable(String infileS, String taxaListFileS, String outfileS){

        List<Integer> indexTaxa = new ArrayList<>();
        String[] taxaArray = AoFile.getStringArraybyList_withoutHeader(taxaListFileS,0);
        Arrays.sort(taxaArray);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中

                }
                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write("CHROM" + "\t" + l.get(1));
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(taxaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                            bw.write("\t" + l.get(i));
                        }
                    }
                    bw.newLine(); //写完之后记得换行
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    String altList = l.get(4);
                    int nAlt = PStringUtils.fastSplit(altList, ",").size();
                    if (nAlt > 1) continue; //filter alt num with 2 or more
                    if (altList.equals("D") || altList.equals("I")) continue; //filter D I

                    String chr = l.get(0);
                    String pos = l.get(1);
                    bw.write(chr + "\t" + pos);

                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        String geno = l.get(indexTaxa.get(i));
                        geno = PStringUtils.fastSplit(geno,":").get(0);
//                        System.out.println(geno);
                        if(geno.equals("0/0")){
                            bw.write( "\t0");
                        }
                        if(geno.equals("0/1") || geno.equals("1/0")){
                            bw.write(  "\t1");
                        }
                        if(geno.equals("1/1")){
                            bw.write( "\t2");
                        }
                        if(geno.equals("./.")){
                            bw.write("\tNA");
                        }
                    }
                    bw.newLine();
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexTaxa.size() + "\tTotal sites : " + cnt);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * 从VCF文件中提取 pos 文件，
     *
     * @param infileS   VCF file
     */
    public static TIntArrayList extractVCFPos(String infileS) {
        TIntArrayList out = new TIntArrayList();
        try {
            BufferedReader br = AoFile.readFile(infileS);

            String temp = null;
            List<String> l = new ArrayList<>();
            int cntRaw = 0; //1.总共的SNP数量

            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("#")) continue;
                //***********************************************************//
                if (!temp.startsWith("#")) {
                    cntRaw++;
                    temp = temp.substring(0,20);
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    out.add(pos);
                } //
            }
            br.close();
//            System.out.println(infileS + " is completed. snpNum is\t" + out.size());
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }

    /**
     * 根据提供的taxa列表，从总的VCF文件中提取所需的 chr pos 文件，并对没有分离的位点进行去除,没有分离位点包括全部都是./.的位点
     *
     * @param infileS   VCF file
     * @param taxalist !!!! 没有header，一行一个taxa名字
     */
    public static TIntArrayList extractVCFchrPos(String infileS, String taxalist) {
        TIntArrayList out = new TIntArrayList();
        List<Integer> indexTaxa = new ArrayList<>();
        String[] taxaArray = AoFile.getStringArraybyList_withoutHeader(taxalist,0);
        System.out.println("Chr\tNum_MergedFileVariants\tNum_KeptVariants\tNum_RemovedSites\tNum_NosegregationSites\tNum_NogenotypeSites");
        try {
            BufferedReader br = AoFile.readFile(infileS);

            String temp = null;
            List<String> l = new ArrayList<>();
            int cntRaw = 0; //1.总共的SNP数量
            int cntKept = 0; //2.提取后保留的SNP数量
            int cntRemoved = 0;  //3.去除的SNP数量
            int cntNosegregation = 0; //4.没有分离位点的sites
            int cntSiteNogeno = 0; //5.没有基因型的sites

            while ((temp = br.readLine()) != null) {
                int cntNogenotype = 0;
                //***********************************************************//
                if (temp.startsWith("##")) continue;
                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(taxaArray, taxon);
                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                        }
                    }
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("#")) {
                    cntRaw++;
                    l = PStringUtils.fastSplit(temp);
                    List<String> lTaxaGeno = new ArrayList<>();
                    String altList = l.get(4);
                    int pos = Integer.parseInt(l.get(1));
                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        lTaxaGeno.add(l.get(indexTaxa.get(i)));
                    }

                    for (int i = 0; i < lTaxaGeno.size(); i++) { //判断没有基因型的taxa数目
                        if (lTaxaGeno.get(i).startsWith(".")) {
                            cntNogenotype++;
                        }
                    }

                    if (cntNogenotype == lTaxaGeno.size()) { //过滤 所有taxa都没有基因型的位点
                        cntSiteNogeno++;
                        continue;
                    } //若不过滤，则全是./.的位点在下面的分离测试中会统计到
                    String[] taxaGenoArray = lTaxaGeno.toArray(new String[lTaxaGeno.size()]);
                    boolean segregation = new CalVCF().ifSegregationIncl2alt(taxaGenoArray, altList);
                    if (segregation == false) { //过滤没有分离的位点
                        cntNosegregation++;
                        continue;
                    }
                    cntKept++;
                    out.add(pos);
                } //
            }
            cntRemoved = cntSiteNogeno + cntNosegregation;
            br.close();
            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntRaw + "\t" + cntKept + "\t" + cntRemoved + "\t" + cntNosegregation + "\t" + cntSiteNogeno);
            System.out.println(infileS + " is completed at " + "\tActual taxa size: " + indexTaxa.size() + "\tGoal taxa size : " + l.size());
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }


    /**
     * 根据提供的taxa列表，从总的VCF文件中提取所需的VCF文件，并对没有分离的位点进行去除,没有分离位点包括全部都是./.的位点
     *
     * @param infileS   VCF file
     * @param outfileS  VCF file
     * @param taxalist !!!! 没有header，一行一个taxa名字
     */
    public void extractVCF(String infileS, String outfileS, String taxalist) {

        List<Integer> indexTaxa = new ArrayList<>();
        String[] taxaArray = AoFile.getStringArraybyList_withoutHeader(taxalist,0);


        System.out.println("Chr\tNum_MergedFileVariants\tNum_KeptVariants\tNum_RemovedSites\tNum_NosegregationSites\tNum_NogenotypeSites");
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cntRaw = 0; //1.总共的SNP数量
            int cntKept = 0; //2.提取后保留的SNP数量
            int cntRemoved = 0;  //3.去除的SNP数量
            int cntNosegregation = 0; //4.没有分离位点的sites
            int cntSiteNogeno = 0; //5.没有基因型的sites

            while ((temp = br.readLine()) != null) {
                int cntNogenotype = 0;
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中
                    bw.write(temp);
                    bw.newLine();
                }
                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(l.get(0));
                    for (int i = 1; i < 9; i++) { //先写前8列的信息
                        bw.write("\t" + l.get(i));
                    }

                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(taxaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexTaxa.add(i);
                            bw.write("\t" + l.get(i));
                        }
                    }
                    bw.newLine(); //写完之后记得换行
                    Collections.sort(indexTaxa);

                }
                if (!temp.startsWith("#")) {
                    cntRaw++;
                    l = PStringUtils.fastSplit(temp);
                    List<String> lTaxaGeno = new ArrayList<>();
                    String altList = l.get(4);
                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        lTaxaGeno.add(l.get(indexTaxa.get(i)));
                    }

                    for (int i = 0; i < lTaxaGeno.size(); i++) { //判断没有基因型的taxa数目
                        if (lTaxaGeno.get(i).startsWith(".")) {
                            cntNogenotype++;
                        }
                    }

                    if (cntNogenotype == lTaxaGeno.size()) { //过滤 所有taxa都没有基因型的位点 全部 ./.
                        cntSiteNogeno++;
                        continue;
                    } //若不过滤，则全是./.的位点在下面的分离测试中会统计到
                    String[] taxaGenoArray = lTaxaGeno.toArray(new String[lTaxaGeno.size()]);
                    boolean segregation = this.ifSegregationIncl2alt(taxaGenoArray, altList);
                    if (segregation == false) { //过滤没有分离的位点
                        cntNosegregation++;
                        continue;
                    }
                    cntKept++;
                    String info = this.getInfo(taxaGenoArray, altList);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 7; i++) {
                        sb.append(l.get(i)).append("\t");
                    }
                    sb.append(info).append("\tGT:AD:PL");
                    for (int i = 0; i < lTaxaGeno.size(); i++) {
                        sb.append("\t").append(lTaxaGeno.get(i));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                } //
            }
            cntRemoved = cntSiteNogeno + cntNosegregation;
            br.close();
            bw.flush();
            bw.close();
            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntRaw + "\t" + cntKept + "\t" + cntRemoved + "\t" + cntNosegregation + "\t" + cntSiteNogeno);
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexTaxa.size() + "\tGoal taxa size : " + l.size());
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false,若全是./.，则返回false,针对有1个或者2个alt的情况，包含D
     * I
     *
     * @param genoArray
     * @param altList
     * @return
     */
    public static boolean ifSegregationIncl2alt(String[] genoArray, String altList) {
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int cntNogenotype = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                cntNogenotype++;
                continue; //表示不再进行下面的计算，跳出该taxa基因型的统计
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
        }
        int cntgeno = genoArray.length - cntNogenotype;
        int cntallele = 2 * cntgeno; //若全部为./.,则cntallele为0,下文也判断为没有分离！！
        //判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false
        boolean a = false;
        if (altList.length() == 1) { //分 alt 是1个的情况
            if ((acCnt[0] == cntallele) || (acCnt[1] == cntallele)) {
                a = false;
            } else {
                a = true;
            }

        } else if (altList.length() > 1) { //alt是2个的情况
            if (acCnt[0] == cntallele || acCnt[1] == cntallele || acCnt[2] == cntallele) {
                a = false;
            } else {
                a = true;
            }
        }
        return a;
    }


    /**
     * 用于获取 geno 的Info信息
     * @param genoArray
     * @param altList
     * @return
     */
    public String getInfo(String[] genoArray, String altList) {
        int dp = 0; //总深度
        int nz = 0; //有基因型的个体数
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的深度统计
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的基因型统计
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合

            //先计算深度
            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
                dp += c; //dp是总深度
                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
            }

            //再计算基因型
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合  如： 0 0
            for (int j = 0; j < temList.size(); j++) { //0/0:13,0:0,4,25
                int c = Integer.parseInt(temList.get(j)); // c是基因型0 1 2 其中的一个
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
            int index1 = Integer.parseInt(temList.get(0)); //0/0基因型的
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
            if (index1 != index2) {
                ht++;
            }
        }
        nz = genoArray.length - nz;
        int sum = 0; //所有allele的总数
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }

        StringBuilder sb = new StringBuilder();

        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1); //删除最后一个字符","号
        sb.append(";AC=");
        for (int i = 1; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) { //二维数组的长度是第一维的长度，这里是2（只有1个alt） 或者3 (有2个alt)
            for (int j = i + 1; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";HT=").append(ht).append(";MAF=").append(String.format("%.4f", maf));
        return sb.toString();
    }


    /**
     * 根据提供的taxa列表，从总的VCF文件中提取所需的VCF文件，并对没有分离的位点进行去除,没有分离位点包括全部都是./.的位点
     *
     * @param infileS
     * @param outfileS
     * @param taxalist 没有header，一行一个taxa名字
     */
    public static void calMafFromPop(String infileS, String outfileS, String taxalist) {
        List<Integer> indexTaxa = new ArrayList<>();
        String[] taxaArray = AoFile.getStringArraybyList_withoutHeader(taxalist,0);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Maf");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(taxaArray, taxon);
                        if (index1 > -1) {
                            indexTaxa.add(i);
                        }
                    }
                    Collections.sort(indexTaxa);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    List<String> lTaxaGeno = new ArrayList<>();
                    for (int i = 0; i < indexTaxa.size(); i++) { //无论有无基因型，都加进去了
                        lTaxaGeno.add(l.get(indexTaxa.get(i)));
                    }
                    String[] taxaGenoArray = lTaxaGeno.toArray(new String[lTaxaGeno.size()]);
                    String maf = CalVCF.getPopMAF(taxaGenoArray);
                    bw.write(maf);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static String getPopMAF(String[] PopGenoArray) {
        int[] acCnt = new int[2]; //所有包括ref和alt的个数
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) continue;
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
        }
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf > 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }

        StringBuilder sb = new StringBuilder();
        if (maf == 0 || maf == 1) {
            sb.append("NA");
        }else{
            sb.append(String.format("%.4f", maf));
        }

        return sb.toString();
    }

    /**
     * extract pos info from vcf file. eg: vcf ---- ID Chr Pos Ref Alt
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public static void extractIDHapPosRefAlt(String infileDirS, String outfileDirS) {

        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".pos.Base.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".pos.Base.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnt = 0;
                bw.write("ID\tChr\tPos\tRef\tAlt\n");
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temp = temp.substring(0, 40); //肯定够
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb = new StringBuilder(l.get(2));
                    sb.append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt % 1000000 == 0) {
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
     * extract pos info from vcf file. eg: vcf ---- Chr Pos
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public static void mkHapPosWithoutHeader(String infileDirS, String outfileDirS) {
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

                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int cnt = 0;
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
                String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + ".pos.txt.gz").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
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
     * 根据barley的fasta格式，写出phylip格式，连续型的
     */
    public void mkPhylipFormat(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/A_subgenome.fasta.txt";
        String outfileS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/chrA_sub_barley_phylip.txt";

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> fasta = new ArrayList<>();
            List<String> l = new ArrayList<>();

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    continue;
                }
                else {
                    for (int i = 0; i < temp.length(); i++) {
                        fasta.add(temp.substring(i,i+1));
                    }
                }
            }
            br.close();

            //确定要写几行
            int lines = Integer.MIN_VALUE;
            if(fasta.size()%50 == 0){
                lines = fasta.size()/50;
            }
            else{
                lines = fasta.size()/50 + 1;
            }

            //先写第一行
            int k = -1;
            bw.write("Barley");
            int space = 15 - "Barley".length();
            for (int i = 0; i < space; i++) {
                bw.write(" ");
            }
            for (int i = 0; i < 50; i++) {
                k++;
                bw.write(fasta.get(k));
            }
            bw.newLine();

            for (int i = 1; i < lines; i++) {
                //开始写第二行
                for (int j = 0; j < 65; j++) {
                    if(j<15){
                        bw.write(" ");
                    }
                    else{
                        k++;
                        if(k<fasta.size()){
                            bw.write(fasta.get(k));
                        }
                    }
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);

        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }


    /**
     * 改变VCF的taxa名字
     *
     */
    public void reheader(String infileS, String outfileS, String reheaderS){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/chr.Asubgenome.vcf.gz";
//        String outfileS= "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/chr.Asubgenome_reheader.vcf.gz";
//        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hm.put(taxa,taxaID);
        }

        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            List<String> taxaList = new ArrayList<>();
            List<String> l = new ArrayList<>();

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    bw.write(temp);
                    bw.newLine();
                }
                else if(temp.startsWith("#C")){
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 9; i++) {
                        sb.append(l.get(i)).append("\t");
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    for (int i = 9; i < l.size(); i++) {
                        String taxa = l.get(i);
                        taxa = hm.get(taxa);
                        sb.append("\t").append(taxa);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                else if(!temp.startsWith("#")) { //
                    bw.write(temp);
                    bw.newLine();

                }//
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);

        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     * 根据提供的taxa列表，从总的VCF文件中提取所需的VCF文件，并对没有分离的位点进行去除,没有分离位点包括全部都是./.的位点
     * 计算特定群体的位点杂合度
     *
     * @param infileS
     * @param outfileS
     * @param taxalist 没有header，一行一个taxa名字
     */
    public void getSNPHeter(String infileS, String outfileS, String taxalist) {
        List<Integer> indexHexa = new ArrayList<>();
        List<String> lhexa = new AoFile().getStringListwithoutHeader(taxalist,0); //六倍体的taxa名集合
        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        Arrays.sort(hexaArray);
        System.out.println("Chr\tNum_MergedFileVariants\tNum_KeptVariants\tNum_RemovedSites\tNum_NosegregationSites\tNum_NogenotypeSites");
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);

            bw.write("Chr\tPos\tHetProportion");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();

            int cntRaw = 0; //1.总共的SNP数量
            int cntKept = 0; //2.提取后保留的SNP数量
            int cntRemoved = 0;  //3.去除的SNP数量
            int cntNosegregation = 0; //4.没有分离位点的sites
            int cntSiteNogeno = 0; //5.没有基因型的sites

            while ((temp = br.readLine()) != null) {
                int cntNogenotype = 0;
                //***********************************************************//
                if (temp.startsWith("##")) continue;//
                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexHexa.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                }
                if (!temp.startsWith("#")) {
                    cntRaw++;
                    l = PStringUtils.fastSplit(temp);
                    List<String> lHexaGeno = new ArrayList<>();
                    String altList = l.get(4);
                    for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                        lHexaGeno.add(l.get(indexHexa.get(i)));
                    }

                    for (int i = 0; i < lHexaGeno.size(); i++) { //判断没有基因型的taxa数目
                        if (lHexaGeno.get(i).startsWith(".")) {
                            cntNogenotype++;
                        }
                    }

                    if (cntNogenotype == lHexaGeno.size()) { //过滤 所有taxa都没有基因型的位点
                        cntSiteNogeno++;
                        continue;
                    } //若不过滤，则全是./.的位点在下面的分离测试中会统计到
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    boolean segregation = this.ifSegregation(hexaGenoArray, altList);
                    if (segregation == false) { //过滤没有分离的位点
                        cntNosegregation++;
                        continue;
                    }
                    cntKept++;
                    double h = this.calSNPSitesHeter(hexaGenoArray);
                    bw.write(l.get(0)+"\t"+l.get(1)+"\t");
                    bw.write(String.format("%.4f",h));
                    bw.newLine();
                } //
            }
            cntRemoved = cntSiteNogeno + cntNosegregation;
            br.close();
            bw.flush();
            bw.close();
            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntRaw + "\t" + cntKept + "\t" + cntRemoved + "\t" + cntNosegregation + "\t" + cntSiteNogeno);
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexHexa.size() + "\tGoal taxa size : " + lhexa.size()+ ".");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false,只针对有1个alt的情况
     *
     * @param genoArray
     * @param altList
     * @return
     */
    public static boolean ifSegregation(String[] genoArray, String altList) {
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
        }
        //判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false
        boolean a = false;
        if ((acCnt[0] == 0 && acCnt[1] > 0) || (acCnt[0] > 0 && acCnt[1] == 0)) {
            a = false;
        } else {
            a = true;
        }
        return a;
    }

    /**
     * 根据MAF过滤VCF
     * @param infileS
     * @param mafRatio filter maf less than ratio
     * @param outfileS
     */

    public static void filterMAFinVCF(String infileS, double mafRatio, String outfileS){
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cntTotal = 0;
            int cntKeep = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                if (!temp.startsWith("#")) {
                    cntTotal++;
                    l = PStringUtils.fastSplit(temp);
                    List<String> lgeno = new ArrayList<>();
                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                    }
                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String maf = CalVCF.getPopMAF(genoArray);
                    if (Double.parseDouble(maf) < mafRatio) continue;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                    cntKeep++;
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("****************************** LOG ******************************");
            System.out.println("Total SNP number    " + cntTotal + "    Kept SNP number " + cntKeep);
            System.out.println(infileS + " is completed at " + outfileS);


        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


    }

    /**
     * 根据位点杂合度过滤单个群体的VCF
     *
     * @param infileS
     * @param outfileS
     */
    public static void filterHeterinVCF(String infileS, double ratio, String outfileS) {

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cntTotal = 0;
            int cntKeep = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                if (!temp.startsWith("#")) {
                    cntTotal++;
                    l = PStringUtils.fastSplit(temp);
                    List<String> lgeno = new ArrayList<>();
                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                    }
                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    double heterRate = CalVCF.calSNPSitesHeter(genoArray);
                    if (heterRate < ratio ) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();
                        cntKeep++;
                    }
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("****************************** LOG ******************************");
            System.out.println("Total SNP number    " + cntTotal + "    Kept SNP number " + cntKeep);
            System.out.println(infileS + " is completed at " + outfileS);


        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * return the site MAF from vcf
     * @param genoArray
     * @return
     */
    public static Double calSNPsitesMAF(String[] genoArray){
        Double out = Double.MIN_VALUE;
        int nz = 0; //有基因型的个体数
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (!genoArray[i].startsWith(".")) {
                nz++;
                tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
                //再计算基因型
                temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合 0/0的集合
                int index1 = Integer.parseInt(temList.get(0)); //0/0基因型的
                int index2 = Integer.parseInt(temList.get(1));
                if (index1 != index2) {
                    ht++;
                }
            }
        }
        out = (double) ht/nz;
        return out;
    }

    /**
     * return the site heterozygosity from vcf
     *
     * @param genoArray
     * @return
     */
    public static Double calSNPSitesHeter(String[] genoArray){
        Double out = Double.MIN_VALUE;
        int nz = 0; //有基因型的个体数
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (!genoArray[i].startsWith(".")) {
                nz++;
                tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
                //再计算基因型
                temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合 0/0的集合
                int index1 = Integer.parseInt(temList.get(0)); //0/0基因型的
                int index2 = Integer.parseInt(temList.get(1));
                if (index1 != index2) {
                    ht++;
                }
            }
        }
        out = (double) ht/nz;
        return out;
    }

//vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --weir-fst-pop China_specific.txt --weir-fst-pop Mixed.txt --out China_specificVSMixed_chr001 &
//vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --weir-fst-pop Non_stiff_stalk.txt --weir-fst-pop China_specific.txt --fst-window-size 10000 --fst-window-step 2000 --out Non_stiff_stalkVSChina_specific_chr001 &
//vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --keep China_specific.txt --site-pi --out China_specific_hmp321_agpv4_chr10 &
//vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --keep China_specific.txt --window-pi 10000 --out China_specific_hmp321_agpv4_chr10 &
//vcftools --site-mean-depth --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/singleChr/chr036.Dlineage.maf0.005.bi_subset.vcf.gz &
// vcftools --depth --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/singleChr/chr036.Dlineage.maf0.005.bi_subset.vcf.gz --out chr036.Dlineage.maf0.005.bi_subset &
    
    public void getSiteMeanDepth() {
        String infileDirS = "";
        String outfileDirS = "";
        List<Integer> lA = new ArrayList<>();
        List<Integer> lD = new ArrayList<>();
        //先进行D的建立
        int j = 5;
        lD.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            lD.add(j);
        }
        int k = 6;
        lD.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            lD.add(k);
        }
        //再进行A的建立
        int a = 1;
        lA.add(a);
        for (int i = 0; i < 6; i++) {
            a = a + 6;
            lA.add(a);
        }
        int aa = 2;
        lA.add(aa);
        for (int i = 0; i < 6; i++) {
            aa = aa + 6;
            lA.add(aa);
        }

        Collections.sort(lA);
        Collections.sort(lD);
        
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(lD, i);
            int index2 = Collections.binarySearch(lA, i);
            if (index < 0) { //说明是属于AB的
                if (index2 > -1) { //说明是属于A的
                    String mPath = new File(infileDirS, "chr" + chr + ".Alineage.vcf").getAbsolutePath();
                    System.out.println("vcftools --site-mean-depth --gzcvf " + mPath );
                } else { //说明是属于B的
                    String mPath = new File(infileDirS, "chr" + chr + ".Blineage.vcf").getAbsolutePath();
                    System.out.println();
                }
            } else if (index > -1) { //说明是属于D的
                String mPath = new File(infileDirS, "chr" + chr + ".Dlineage.vcf").getAbsolutePath();
                System.out.println();
            }
        }
    }
    
    public void bcftools_merge() {
        String abdFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/output/VCF/";
        String abFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/ab/output/VCF/";
        String dFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/d/output/VCF/";
        String mergedFileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/";
        /**
         * pseudo-code: 1.建立3个lineage的list,然后进行循环，判断：在A lineage下合并，依次类推。
         */
        List<Integer> lA = new ArrayList<>();
        List<Integer> lD = new ArrayList<>();
        //先进行D的建立
        int j = 5;
        lD.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            lD.add(j);
        }
        int k = 6;
        lD.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            lD.add(k);
        }
        //再进行A的建立
        int a = 1;
        lA.add(a);
        for (int i = 0; i < 6; i++) {
            a = a + 6;
            lA.add(a);
        }
        int aa = 2;
        lA.add(aa);
        for (int i = 0; i < 6; i++) {
            aa = aa + 6;
            lA.add(aa);
        }

        Collections.sort(lA);
        Collections.sort(lD);
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            String abdPath = new File(abdFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String abPath = new File(abFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String dPath = new File(dFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            int index = Collections.binarySearch(lD, i);
            int index2 = Collections.binarySearch(lA, i);
            if (index < 0) { //说明是属于AB的
                if (index2 > -1) { //说明是属于A的
                    String mPath = new File(mergedFileDirS, "chr" + chr + ".Alineage.vcf").getAbsolutePath();
                    System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + abPath + " -o " + mPath + " &");
                } else { //说明是属于B的
                    String mPath = new File(mergedFileDirS, "chr" + chr + ".Blineage.vcf").getAbsolutePath();
                    System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + abPath + " -o " + mPath + " &");
                }
            } else if (index > -1) { //说明是属于D的
                String mPath = new File(mergedFileDirS, "chr" + chr + ".Dlineage.vcf").getAbsolutePath();
                System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + dPath + " -o " + mPath + " &");
            }
        }
    }

    public void bgzip_ABD() {
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
        }
    }

}
