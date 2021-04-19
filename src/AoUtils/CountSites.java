/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class CountSites {

    public CountSites() {
//        this.getSharedSNP();
        //this.mergesubsetVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/Asub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/004_merged/chr.Asubgenome.maf0.01.SNP_bi.subset.vcf.gz");
        //this.mergesubsetVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/Bsub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/004_merged/chr.Bsubgenome.maf0.01.SNP_bi.subset.vcf.gz");
        //this.mergesubsetVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/Dsub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/004_merged/chr.Dsubgenome.maf0.01.SNP_bi.subset.vcf.gz");
        //this.mergesubsetVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/004_merged","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/005_all/chr.ABsubgenome.maf0.01.SNP_bi.subset.vcf.gz");
        //this.calSNPHetMissMaf("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/001_rawVCF", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/002_calMAF");
        //new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/006_calMAF", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/007_bintable", "25", "0.5");
        //this.mergesubsetVCF("", "");
//        this.extractVCF("/Users/Aoyue/Documents/chr001.ABDgenome.filtered0.75_10000lines.vcf", "/Users/Aoyue/Documents/test.vcf", "/Users/Aoyue/Documents/a.txt");

//        this.cntSitesinMergedVCFtoPop("/Users/Aoyue/Documents/a/chr036.vcf", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt");
//        this.mergeChr1and2txt_double("/Users/Aoyue/Documents/a.txt", "/Users/Aoyue/Documents/b.txt");
//        this.getRefChrPos(2, 2);
//        mergefileandChangeChrPos_chr1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/001_ori/delSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos");


//        this.changechrPosOnVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/008_VCF/001_exonVCF/000_exonVCF");
//    this.mergeVCFfile1and2_chr1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/008_VCF/001_exonVCF/000_exonVCF","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/008_VCF/001_exonVCF/001_mergeVCF");
//        this.mergeChr1and2txt_int("/Users/Aoyue/Downloads/log_abd_countSitesinFastCallformat_20200519.txt","/Users/Aoyue/Downloads/log_abd.txt");
//        this.mergeChr1and2txt_int("/Users/Aoyue/Downloads/log_maf0.01_countSitesinFastCallformat20200520.txt","/Users/Aoyue/Downloads/log_maf0.01_d_countSitesinFastCallformat20200520.txt");
//        this.mergeChr1and2_Dgenome("/Users/Aoyue/Downloads/log_maf0.01_countSitesinFastCallformat20200520.txt","/Users/Aoyue/Downloads/log_maf0.01_ddd_countSitesinFastCallformat20200520.txt");


    }

    /**
     * 从Tassel软件的distance matrix结果中提取某一行的结果，即所有taxa到指定ref taxa的遗传距离
     * ##IBS_Distance_Matrix.AverageTotalSites=444322.6751757548
     * ##IBS_Distance_Matrix.NumAlleles=3
     * ##IBS_Distance_Matrix.TrueIBS=false
     * ##Matrix_Type=IBS_Distance_Matrix
     * 607
     */
    public static void extractIBSdistanceFromMatrix(String infileS, String outfileS,String Reftaxa){
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //过滤前5行，并打印taxa数
            String temp = null;
            int cntTaxaNum = 0;
            for (int i = 0; i < 5; i++) {
                temp = br.readLine();
                if (i==4) {
                    cntTaxaNum = Integer.parseInt(temp);
                    System.out.println("Total taxa num\t" + cntTaxaNum);
                }
            }
            List<String> l = new ArrayList<>();
            List<String> taxaList = new ArrayList<>();
            List<String> valueList = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(0);
                taxaList.add(l.get(0));
                if (taxa.equals(Reftaxa)){
                    for (int i = 1; i < l.size(); i++) {
                        valueList.add(l.get(i));
                    }
                }
                cnt++;
            }
            System.out.println("Total taxaList num\t" + taxaList.size());
            br.close();

            bw.write("Taxa\tIBSdistance");bw.newLine();
            for (int i = 0; i < taxaList.size(); i++) {
                bw.write(taxaList.get(i) + "\t" + valueList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }



    /**
     * 根据合并后的VCF文件，进行群体内部的计数，为了得到群体的SNP和Indel以及未分离的位点数目。
     *
     * @param infileS
     * @param taxalist
     */
    public void cntSitesinMergedVCFtoPop(String infileS, String taxalist) {
        List<String> lhexa = new ArrayList<>(); //六倍体的taxa名集合
        List<Integer> indexHexa = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(taxalist);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        Arrays.sort(hexaArray);
        System.out.println("Chr\tRaw Variants\tKept Variants\tNogenotype Sites\tNosegregation Sites\tTotal SNPs\tBiallelic SNPs\tTriallelic SNPs\tIndels\tInsertions\tDeletions");
//        System.out.println("Chr\tRaw Variants\tKept Variants\tNosegregation Sites");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
//            if (outfileS.endsWith(".vcf")) {
//                bw = IOUtils.getTextWriter(outfileS);
//            } else if (outfileS.endsWith(".vcf.gz")) {
//                bw = IOUtils.getTextGzipWriter(outfileS);
//            }

            String temp = null;
            List<String> l = new ArrayList<>();
            String count = null;

            int cntRaw = 0; //
            int cntKept = 0; //
            int cntSiteNogeno = 0; //
            int cntNosegregation = 0;

            int cntSNP = 0;
            int cntBi = 0;
            int cntTri = 0;
            int cntIndel = 0;
            int cntI = 0;
            int cntD = 0;

            while ((temp = br.readLine()) != null) {
                int cntNogenotype = 0;
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中
//                    bw.write(temp);
//                    bw.newLine();
                }
                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
//                    bw.write(l.get(0));
//                    for (int i = 1; i < 9; i++) { //先写前8列的信息
//                        bw.write("\t" + l.get(i));
//                    }

                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);

                        if (index1 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexHexa.add(i);
//                            bw.write("\t" + l.get(i));
                        }
                    }
//                    bw.newLine(); //写完之后记得换行
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
                    } //若不过滤，则全是./.的位点在下面的分离测试中会加一
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    boolean segregation = Boolean.valueOf(this.ifSegregationandrealAlt(hexaGenoArray, altList).split(";")[0]);
                    if (segregation == false) { //过滤没有分离的位点
                        cntNosegregation++;
                        continue;
                    }
                    cntKept++;

                    String alt = this.ifSegregationandrealAlt(hexaGenoArray, altList).split(";")[1];
                    if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                        boolean ifD = false;
                        if (!alt.contains("D") && (!alt.contains("I"))) {
                            cntTri++;
                            cntSNP++;
                        }
                        if (alt.contains("D")) {
                            cntD++;
                            cntIndel++;
                            ifD = true;
                        }
                        if (alt.contains("I")) {
                            cntI++;
                            if (ifD == false) {
                                cntIndel++;
                            }
                        }

                    } else if (alt.length() == 1) { //1个alt的情况;
                        if (!alt.equals("D") && (!alt.equals("I"))) {
                            cntBi++;
                            cntSNP++;
                        }
                        if (alt.equals("D")) {
                            cntD++;
                            cntIndel++;
                        }
                        if (alt.equals("I")) {
                            cntI++;
                            cntIndel++;
                        }
                    }

//                    StringBuilder sb = new StringBuilder();
//                    for (int i = 0; i < 7; i++) {
//                        sb.append(l.get(i)).append("\t");
//                    }
//                    sb.append(".").append("\tGT:AD:PL");
//                    for (int i = 0; i < lHexaGeno.size(); i++) {
//                        sb.append("\t").append(lHexaGeno.get(i));
//                    }
//                    bw.write(sb.toString());
//                    bw.newLine();
                } //
            }
//            br.close();
//            bw.flush();
//            bw.close();

            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntRaw + "\t" + cntKept + "\t" + cntSiteNogeno + "\t" + cntNosegregation
                    + "\t" + cntSNP + "\t" + cntBi + "\t" + cntTri + "\t" + cntIndel + "\t" + cntI + "\t" + cntD);
            System.out.println(infileS + " is completed. " + "\tActual taxasize: " + indexHexa.size() + "\tGoal taxa size : " + lhexa.size());

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false,若全是./.，则返回false,针对有1个或者2个alt的情况，包含D
     * I,并且返回所抽取样品的ALT信息
     *
     * @param genoArray
     * @param altList
     * @return
     */
    public String ifSegregationandrealAlt(String[] genoArray, String altList) {
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
        String alt = altList;
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
                if (acCnt[1] == 0) //第一个alt是空的
                {
                    alt = altList.split(",")[1]; //返回第二个alt
                }
            }
            if (acCnt[2] == 0) { //第2个alt空的
                alt = altList.split(",")[0];
            }
        }

        return String.valueOf(a) + ";" + alt;
    }


    public void getSharedSNP() {
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr005.ABDgenome.10000lines.vcf";
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr005.Dgenome.10000lines.vcf";
//        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr005.ABDgenome.10000lines.shared.txt";
//        String outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr005.Dgenome.10000lines.shared.txt";
//        String infileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/test.vcf";
//        String outfileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/test.all.shared.txt";

//////// 1A merge test
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr001.ABDgenome.10000lines.vcf";
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr001.ABgenome.10000lines.vcf";
//        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr001.ABDgenome.10000lines.shared.txt";
//        String outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr001.ABgenome.10000lines.shared.txt";
//        String infileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/chr001_mergeTaxa.vcf";
//        String outfileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr001_mergeTaxa.all.shared.txt";
/////// 1B merge test
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr003.ABDgenome.10000lines.vcf";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr003.ABgenome.10000lines.vcf";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr003.ABDgenome.10000lines.shared.txt";
        String outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr003.ABgenome.10000lines.shared.txt";
        String infileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/chr003_mergeTaxa.vcf";
        String outfileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/chr003_mergeTaxa.all.shared.txt";

        List<String> posl = new ArrayList<>();
        List<String> sharedpos = new ArrayList<>();
        try {
            //先建立pos数据库
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedWriter bw = null;
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                String chr = l.get(0);
                String pos = l.get(1);
                String ref = l.get(3);
                String alt = l.get(4);
                posl.add(pos);
                cnt++;
            }
            br.close();
            System.out.println(cnt + "  snps totally in chr005ABD");

            //找出D中共有的snp,并写出
            Collections.sort(posl);
            br = IOUtils.getTextReader(infileS2);
            bw = IOUtils.getTextWriter(outfileS2);
            bw.write("CHROM\tPOS\tREF\tALT\n");
            int share = 0;
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                String chr = l.get(0);
                String pos = l.get(1);
                String ref = l.get(3);
                String alt = l.get(4);
                int index = Collections.binarySearch(posl, pos);
                if (index >= 0) {
                    share++;
                    sharedpos.add(pos);
                    bw.write(chr + "\t" + pos + "\t" + ref + "\t" + alt);
                    bw.newLine();
                }
            }

            bw.flush();
            bw.close();
            br.close();
            System.out.println(share + "  shared snps totally in chr005ABD and chr005D");

            //写出ABD中共有的snp
            Collections.sort(sharedpos);
            br = IOUtils.getTextReader(infileS1);
            bw = IOUtils.getTextWriter(outfileS1);
            bw.write("CHROM\tPOS\tREF\tALT\n");
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                String chr = l.get(0);
                String pos = l.get(1);
                String ref = l.get(3);
                String alt = l.get(4);
                int index = Collections.binarySearch(sharedpos, pos);
                if (index >= 0) {
                    bw.write(chr + "\t" + pos + "\t" + ref + "\t" + alt);
                    bw.newLine();
                }
            }

            bw.flush();
            bw.close();
            br.close();

            //找出合并后的文件中，共有POS的ALT变化情况
            Collections.sort(sharedpos);
            br = IOUtils.getTextReader(infileS3);
            bw = IOUtils.getTextWriter(outfileS3);
            bw.write("CHROM\tPOS\tREF\tALT\n");
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                String chr = l.get(0);
                String pos = l.get(1);
                String ref = l.get(3);
                String alt = l.get(4);
                int index = Collections.binarySearch(sharedpos, pos);
                if (index >= 0) {
                    bw.write(chr + "\t" + pos + "\t" + ref + "\t" + alt);
                    bw.newLine();
                }
            }

            bw.flush();
            bw.close();
            br.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 将抽样的28条染色体按照Asub Bsub 合并
     * 必须是28条染色体都存在
     * @param infileDirS
     * @param outfileDirS
     */
    public void mergeVCFtoAandBsubgenome(String infileDirS, String outfileDirS) {

        //建立1-42 一一对应的亚基因组的关系，根据chr001找到chr.Asub
        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        HashMap<Integer, Integer> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i],0);
            hml.put(arrb[i], 1);
        }

        //列出文件,并加入对应的文件列表 ## 解释：3个数组，一个数组是一个文件列表，里面含有每个亚基因组对应的文件列表
        List<File>[] fs = new List[2];
        for (int i = 0; i < fs.length; i++) { // 一定记得初始化
            fs[i] = new ArrayList<>();
        }
        File[] fall = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fall);
        for (int i = 0; i < fall.length; i++) {
            String name = fall[i].getName().substring(3,6);
            int chr = Integer.parseInt(name);
            int index = hml.get(chr);
            fs[index].add(fall[i]);
        }

        //设定输出文件的路径及名字
        String[] outfileS = new String[2];
        outfileS[0] = new File(outfileDirS, "Asubgenome.vcf").getAbsolutePath();
        outfileS[1] = new File(outfileDirS, "Bsubgenome.vcf").getAbsolutePath();

        //开始进行写文件
        BufferedWriter[] bws = new BufferedWriter[2];
        for (int i = 0; i < bws.length; i++) {
            bws[i] = AoFile.writeFile(outfileS[i]);
        }

        //根据亚基因组文件列表，确定输出文件的index
        HashMap<List<File>,BufferedWriter> hmFilesubOutS = new HashMap<>();
        hmFilesubOutS.put(fs[0],bws[0]);
        hmFilesubOutS.put(fs[1],bws[1]);


        List<List<File>> fsSubList = Arrays.asList(fs); //表明数组是list类型，list又是file类型
        fsSubList.stream().forEach(p ->{ // p是一个list

            try{
                //先写表头
                String infileS = p.get(0).getAbsolutePath();
                BufferedWriter bw = hmFilesubOutS.get(p);
                BufferedReader br = AoFile.readFile(infileS);
                String temp;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                br.close();

                int total = 0;
                for (int i = 0; i < p.size(); i++) {
                    br = AoFile.readFile(p.get(i).getAbsolutePath());
                    int cnt =0;
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
                    System.out.println(cnt + "\tsnps in " + p.get(i).getName());

                }
                System.out.println(total + "\tsnps totally, mergevcf pipeline is completed\t");
                System.out.println("************************************************************");
                br.close();
                bw.flush();
                bw.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);

            }

        });
    }


    /**
     *
     *  将抽样的42条染色体按照ABsub合并,忽略 5 6 11 12 等D基因组的文件，
     *  必须是42条染色体都存在!!!
     * @param infileDirS
     * @param outfileDirS
     */
    public static void mergeVCFtoABsubgenome(String infileDirS, String outfileDirS) {

        //建立1-42 一一对应的亚基因组的关系，根据chr001找到chr.Asub
        int[] arrab = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38,3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, Integer> hml = new HashMap<>();
        Arrays.sort(arrab);
        Arrays.sort(arrd);
        for (int i = 0; i < arrab.length; i++) {
            hml.put(arrab[i],0);
        }
//        for (int i = 0; i < arrd.length; i++) {
//            hml.put(arrd[i], 1);

//        }

        //列出文件,并加入对应的文件列表 ## 解释：3个数组，一个数组是一个文件列表，里面含有每个亚基因组对应的文件列表
        List<File>[] fs = new List[1];
        for (int i = 0; i < fs.length; i++) { // 一定记得初始化
            fs[i] = new ArrayList<>();
        }
        File[] fall = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fall);
        for (int i = 0; i < fall.length; i++) {
            String name = fall[i].getName().substring(3,6);
            int chr = Integer.parseInt(name);
            /**
             * 去掉D基因组14个文件的合并
             */
            int indexd = Arrays.binarySearch(arrd,chr);
            if (indexd > -1)continue;
            int index = hml.get(chr);
            fs[index].add(fall[i]);
        }

        //设定输出文件的路径及名字
        String[] outfileS = new String[1];
        outfileS[0] = new File(outfileDirS, "ABsubgenome.vcf").getAbsolutePath();
//        outfileS[1] = new File(outfileDirS, "Dsubgenome.vcf").getAbsolutePath();

        //开始进行写文件
        BufferedWriter[] bws = new BufferedWriter[1];
        for (int i = 0; i < bws.length; i++) {
            bws[i] = AoFile.writeFile(outfileS[i]);
        }

        //根据亚基因组文件列表，确定输出文件的index
        HashMap<List<File>,BufferedWriter> hmFilesubOutS = new HashMap<>();
        hmFilesubOutS.put(fs[0],bws[0]);
//        hmFilesubOutS.put(fs[1],bws[1]);


        List<List<File>> fsSubList = Arrays.asList(fs); //表明数组是list类型，list又是file类型
        fsSubList.stream().forEach(p ->{ // p是一个list

            try{
                //先写表头
                String infileS = p.get(0).getAbsolutePath();
                BufferedWriter bw = hmFilesubOutS.get(p);
                BufferedReader br = AoFile.readFile(infileS);
                String temp;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                br.close();

                int total = 0;
                for (int i = 0; i < p.size(); i++) {
                    br = AoFile.readFile(p.get(i).getAbsolutePath());
                    int cnt =0;
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
                    System.out.println(cnt + "\tsnps in " + p.get(i).getName());

                }
                System.out.println(total + "\tsnps totally, mergevcf pipeline is completed\t");
                System.out.println("************************************************************");
                br.close();
                bw.flush();
                bw.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);

            }

        });
    }


    /**
     *
     *  将抽样的42条染色体按照ABsub Dsub 合并
     *  必须是42条染色体都存在
     * @param infileDirS
     * @param outfileDirS
     */
    public static void mergeVCFtoAB_Dsubgenome(String infileDirS, String outfileDirS) {

        //建立1-42 一一对应的亚基因组的关系，根据chr001找到chr.Asub
        int[] arrab = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38,3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, Integer> hml = new HashMap<>();
        Arrays.sort(arrab);
        Arrays.sort(arrd);
        for (int i = 0; i < arrab.length; i++) {
            hml.put(arrab[i],0);
        }
        for (int i = 0; i < arrd.length; i++) {
            hml.put(arrd[i], 1);

        }

        //列出文件,并加入对应的文件列表 ## 解释：3个数组，一个数组是一个文件列表，里面含有每个亚基因组对应的文件列表
        List<File>[] fs = new List[2];
        for (int i = 0; i < fs.length; i++) { // 一定记得初始化
            fs[i] = new ArrayList<>();
        }
        File[] fall = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fall);
        for (int i = 0; i < fall.length; i++) {
            String name = fall[i].getName().substring(3,6);
            int chr = Integer.parseInt(name);
            int index = hml.get(chr);
            fs[index].add(fall[i]);
        }

        //设定输出文件的路径及名字
        String[] outfileS = new String[2];
        outfileS[0] = new File(outfileDirS, "ABsubgenome.vcf").getAbsolutePath();
        outfileS[1] = new File(outfileDirS, "Dsubgenome.vcf").getAbsolutePath();

        //开始进行写文件
        BufferedWriter[] bws = new BufferedWriter[2];
        for (int i = 0; i < bws.length; i++) {
            bws[i] = AoFile.writeFile(outfileS[i]);
        }

        //根据亚基因组文件列表，确定输出文件的index
        HashMap<List<File>,BufferedWriter> hmFilesubOutS = new HashMap<>();
        hmFilesubOutS.put(fs[0],bws[0]);
        hmFilesubOutS.put(fs[1],bws[1]);


        List<List<File>> fsSubList = Arrays.asList(fs); //表明数组是list类型，list又是file类型
        fsSubList.stream().forEach(p ->{ // p是一个list

            try{
                //先写表头
                String infileS = p.get(0).getAbsolutePath();
                BufferedWriter bw = hmFilesubOutS.get(p);
                BufferedReader br = AoFile.readFile(infileS);
                String temp;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                br.close();

                int total = 0;
                for (int i = 0; i < p.size(); i++) {
                    br = AoFile.readFile(p.get(i).getAbsolutePath());
                    int cnt =0;
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
                    System.out.println(cnt + "\tsnps in " + p.get(i).getName());

                }
                System.out.println(total + "\tsnps totally, mergevcf pipeline is completed\t");
                System.out.println("************************************************************");
                br.close();
                bw.flush();
                bw.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);

            }

        });
    }

    /**
     * 将抽样的42条染色体按照Asub Bsub Dsub 合并
     * 必须是42条染色体都存在
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public static void mergeVCFbysubgenome(String infileDirS, String outfileDirS) {

        //建立1-42 一一对应的亚基因组的关系，根据chr001找到chr.Asub
        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, Integer> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i],0);
            hml.put(arrb[i], 1);
            hml.put(arrd[i], 2);
        }

        //列出文件,并加入对应的文件列表 ## 解释：3个数组，一个数组是一个文件列表，里面含有每个亚基因组对应的文件列表
        List<File>[] fs = new List[3];
        for (int i = 0; i < fs.length; i++) { // 一定记得初始化
            fs[i] = new ArrayList<>();
        }
        File[] fall = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fall);
        for (int i = 0; i < fall.length; i++) {
            String name = fall[i].getName().substring(3,6);
            int chr = Integer.parseInt(name);
            int index = hml.get(chr);
            fs[index].add(fall[i]);
        }

        //设定输出文件的路径及名字
        String[] outfileS = new String[3];
        outfileS[0] = new File(outfileDirS, "Asubgenome.vcf").getAbsolutePath();
        outfileS[1] = new File(outfileDirS, "Bsubgenome.vcf").getAbsolutePath();
        outfileS[2] = new File(outfileDirS, "Dsubgenome.vcf").getAbsolutePath();

        //开始进行写文件
        BufferedWriter[] bws = new BufferedWriter[3];
        for (int i = 0; i < bws.length; i++) {
            bws[i] = AoFile.writeFile(outfileS[i]);
        }

        //根据亚基因组文件列表，确定输出文件的index
        HashMap<List<File>,BufferedWriter> hmFilesubOutS = new HashMap<>();
        hmFilesubOutS.put(fs[0],bws[0]);
        hmFilesubOutS.put(fs[1],bws[1]);
        hmFilesubOutS.put(fs[2],bws[2]);


        List<List<File>> fsSubList = Arrays.asList(fs); //表明数组是list类型，list又是file类型
        fsSubList.stream().forEach(p ->{ // p是一个list

            try{
                //先写表头
                String infileS = p.get(0).getAbsolutePath();
                BufferedWriter bw = hmFilesubOutS.get(p);
                BufferedReader br = AoFile.readFile(infileS);
                String temp;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                br.close();

                int total = 0;
                for (int i = 0; i < p.size(); i++) {
                    br = AoFile.readFile(p.get(i).getAbsolutePath());
                    int cnt =0;
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
                    System.out.println(cnt + "\tsnps in " + p.get(i).getName());

                }
                System.out.println(total + "\tsnps totally, mergevcf pipeline is completed\t");
                System.out.println("************************************************************");
                br.close();
                bw.flush();
                bw.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);

            }

        });
    }


    /**
     * 1.将改变位置的VCF文件 chr5和chr6文本文件合并成1D一个文件， chr11,chr12两个文件合并成2D一个文件。自动识别染色体序号并进行合并。
     * 2.针对VCF文件，将其坐标改为参考基因组的坐标 即1A模式
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mergeVCFfileandChangeChrPos_chr1and2(String infileDirS, String outfileDirS) {

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

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
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
                            BufferedReader br = AoFile.readFile(infileS);
                            String outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + fs[i].getName().substring(6).split(".vcf")[0] + ".vcf.gz").getAbsolutePath();
                            //确定输出文件的路径，并读入header
                            String secondchr = PStringUtils.getNDigitNumber(3, chr + 1);
                            //名字变一下：

                            BufferedWriter bw = AoFile.writeFile(outfileS);

                            //读注释信息部分
                            String temp = null;
                            List<String> l = new ArrayList<>();
                            int pos = 0;
                            int a=0;
                            while ((temp = br.readLine()) != null) {
                                if (temp.startsWith("#")) {
                                    bw.write(temp);
                                    bw.newLine();
                                }
                                else{
                                    l = PStringUtils.fastSplit(temp);
                                    StringBuilder sb = new StringBuilder();
                                    chr = Integer.valueOf(l.get(0));
                                    pos = Integer.valueOf(l.get(1));
                                    String Chr = RefV1Utils.getChromosome(chr,pos);
                                    if (chr == 43 || chr == 44) {

                                    } else {
                                        if (chr % 2 == 0) {
                                            pos = RefV1Utils.getPosOnChromosome(chr,pos);
                                        }
                                        sb = new StringBuilder();
                                        sb.append(Chr).append("\t").append(String.valueOf(pos));
                                        for (int k = 2; k < l.size(); k++) {
                                            sb.append("\t");
                                            sb.append(l.get(k));
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

                            ///开始合并文件1和2

                            br.close();

                            //开始读入第2个文件
                            infileS = new File(fs[i].getParent(), fs[i].getName().replaceFirst(PStringUtils.getNDigitNumber(3, chr), secondchr)).getAbsolutePath();
                            br = AoFile.readFile(infileS);
                            while ((temp = br.readLine()) != null) {
                                if (temp.startsWith("#")) {
                                    continue;
                                }
                                l = PStringUtils.fastSplit(temp);
                                StringBuilder sb = new StringBuilder();
                                chr = Integer.valueOf(l.get(0));
                                pos = Integer.valueOf(l.get(1));
                                String Chr = RefV1Utils.getChromosome(chr,pos);
                                if (chr == 43 || chr == 44) {

                                } else {
                                    if (chr % 2 == 0) {
                                        pos = RefV1Utils.getPosOnChromosome(chr,pos);
                                    }
                                    sb = new StringBuilder();
                                    sb.append(Chr).append("\t").append(String.valueOf(pos));
                                    for (int k = 2; k < l.size(); k++) {
                                        sb.append("\t");
                                        sb.append(l.get(k));
                                    }
                                    bw.write(sb.toString());
                                    bw.newLine();
                                    if (a % 1000000 == 0) {
                                        System.out.println("Output " + String.valueOf(a) + " SNPs");
                                    }
                                    a++;
                                }

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
     * 将改变位置的VCF文件 chr5和chr6文本文件合并成1D一个文件， chr11,chr12两个文件合并成2D一个文件。自动识别染色体序号并进行合并。
     * 针对Txt文件
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mergeVCFfile1and2_chr1and2(String infileDirS, String outfileDirS) {

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
                            String outfileS = null;
                            if (infileS.endsWith(".vcf")) {
                                br = IOUtils.getTextReader(infileS);
                                outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + fs[i].getName().substring(6) + ".gz").getAbsolutePath();
                            } else if (infileS.endsWith(".vcf.gz")) {
                                br = IOUtils.getTextGzipReader(infileS);
                                outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + fs[i].getName().substring(6)).getAbsolutePath();
                            }

                            //确定输出文件的路径，并读入header
                            String secondchr = PStringUtils.getNDigitNumber(3, chr + 1);
                            //名字变一下：

                            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);

                            //读注释信息部分
                            String temp = null;
                            while ((temp = br.readLine()) != null) {
                                bw.write(temp);
                                bw.newLine();
                            }

                            ///开始合并文件1和2

                            br.close();

                            //开始读入第2个文件
                            infileS = new File(fs[i].getParent(), fs[i].getName().replaceFirst(PStringUtils.getNDigitNumber(3, chr), secondchr)).getAbsolutePath();
                            if (infileS.endsWith(".vcf")) {
                                br = IOUtils.getTextReader(infileS);
                            } else if (infileS.endsWith(".vcf.gz")) {
                                br = IOUtils.getTextGzipReader(infileS);
                            }

//                            temp = br.readLine(); //read header
                            while ((temp = br.readLine()) != null) {
                                if (temp.startsWith("#")) {
                                    continue;
                                }
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
     * 将改变位置的chr5和chr6文本文件合并成1D一个文件， chr11,chr12两个文件合并成2D一个文件。自动识别染色体序号并进行合并。
     * 针对Txt文件,并且根据Ref改变其CHR POS
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mergefileandChangeChrPos_chr1and2(String infileDirS, String outfileDirS) {

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

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
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
                            BufferedReader br = AoFile.readFile(infileS);
                            String outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + fs[i].getName().substring(6)).getAbsolutePath();

                            //确定输出文件的路径，并读入header
                            String secondchr = PStringUtils.getNDigitNumber(3, chr + 1);
                            //名字变一下：

                            BufferedWriter bw = AoFile.writeFile(outfileS);
                            bw.write(br.readLine()); //先读表头
                            bw.newLine();
                            ///开始合并文件1和2

                            String temp = null; //read header
                            List<String> l = new ArrayList<>();
                            while ((temp = br.readLine()) != null) {
                                l = PStringUtils.fastSplit(temp);
                                String CHR = l.get(0);
                                String POS = l.get(1);
                                String[] change = this.getRefChrPos(Integer.parseInt(CHR), Integer.parseInt(POS));
                                StringBuilder sb = new StringBuilder();
                                sb.append(change[0]).append("\t").append(change[1]);
                                for(int k = 2; k < l.size(); k++){
                                    sb.append("\t").append(l.get(k));
                                }
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                            int a = 3;
                            //开始读入第2个文件
                            infileS = new File(fs[i].getParent(), fs[i].getName().replaceFirst(PStringUtils.getNDigitNumber(3, chr), secondchr)).getAbsolutePath();
                            br = AoFile.readFile(infileS);

                            temp = br.readLine(); //read header
                            while ((temp = br.readLine()) != null) {
                                l = PStringUtils.fastSplit(temp);
                                String CHR = l.get(0);
                                String POS = l.get(1);
                                String[] change = this.getRefChrPos(Integer.parseInt(CHR), Integer.parseInt(POS));
                                StringBuilder sb = new StringBuilder();
                                sb.append(change[0]).append("\t").append(change[1]);
                                for(int k = 2; k < l.size(); k++){
                                    sb.append("\t").append(l.get(k));
                                }
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                            br.close();
                            bw.flush();
                            bw.close();
                            System.out.println(fs[i] + " is completed at " + outfileS);

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
     * 将改变位置的chr5和chr6文本文件合并成1D一个文件， chr11,chr12两个文件合并成2D一个文件。自动识别染色体序号并进行合并。
     * 针对Txt文件
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mergefile1and2_chr1and2(String infileDirS, String outfileDirS) {

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
                            String outfileS = null;
                            if (infileS.endsWith(".txt")) {
                                br = IOUtils.getTextReader(infileS);
                                outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + fs[i].getName().substring(6) + ".gz").getAbsolutePath();
                            } else if (infileS.endsWith(".txt.gz")) {
                                br = IOUtils.getTextGzipReader(infileS);
                                outfileS = new File(outfileDirS, "chr" + hmcntchr.get(chr) + fs[i].getName().substring(6)).getAbsolutePath();
                            }

                            //确定输出文件的路径，并读入header
                            String secondchr = PStringUtils.getNDigitNumber(3, chr + 1);
                            //名字变一下：

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
                            int a = 3;
                            //开始读入第2个文件
                            infileS = new File(fs[i].getParent(), fs[i].getName().replaceFirst(PStringUtils.getNDigitNumber(3, chr), secondchr)).getAbsolutePath();
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
    public void mergeTxtbysuffix(String infileDirS, String outfileS, String suffix) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
//        fs = IOUtils.listFilesEndsWith(fs, suffix);
        fs = IOUtils.listFilesContains(fs, suffix);
        Arrays.sort(fs);
        //System.out.println("Chr\tSNP_Num");

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            ///读表头
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            bw.write(br.readLine());
            bw.newLine();
            
            
            int cnttotal = 0;
            //读正文部分
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
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 目的：将所有txt文本的chr pos位点等信息合并成一个文件。
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
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else if (infileS.endsWith(".fst")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".pi")) {
                br = IOUtils.getTextReader(infileS);
            }

            ///读表头
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            bw.write(br.readLine());
            bw.newLine();

            int cnttotal = 0;
            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }else if (infileS.endsWith(".fst")) {
                    br = IOUtils.getTextReader(infileS);
                }else if (infileS.endsWith(".pi")) {
                    br = IOUtils.getTextReader(infileS);
                }
                String temp = br.readLine(); //read header
//                String temp = null; //not read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
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
     *
     * @param chr
     * @param pos
     * @return
     */
     public String[] getRefChrPos(int chr, int pos) {
        String[] out = new String[2];

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
        String Chr = hmcntchr.get(chr);
        if (chr % 2 == 0) {
            pos = pos + hmcntsPOS.get(chr);
        }
        out[0] = Chr;
        out[1] = String.valueOf(pos);
//        System.out.println(out[0] + " " + out[1]);
        return out;
    }

    /**
     * Change the txt file of chr pos column into 1A mode. Chr Pos 6	291 ---> 1D
     * 452179895
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void changechrPosonTxt(String infileDirS, String outfileDirS) {
        

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
                if (!(i % 2 == 0)) { //i为单数的时候,染色体号为双数
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
            String temp = null;
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
                
                temp = br.readLine();
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
                        for (int i = 2; i < l.size(); i++) {
                            sb.append("\t" + l.get(i));
                        }
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
                System.out.println(temp);
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
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
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
                    temp = temp.substring(0, 50); //肯定够                 
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
        fsList.parallelStream().forEach(f -> {
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
                    if (cnt % 100000 == 0) {
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
     * 结果只有一列 depth平均深度
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void calVcfAverageDepth(String infileDirS, String outfileDirS) {
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
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".depth.txt.gz")).getAbsolutePath();
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".depth.txt.gz")).getAbsolutePath();
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> linetaxa = PStringUtils.fastSplit(temp, "\t");
                bw.write("AverageDepth");
                bw.newLine();
                String[] taxa = new String[linetaxa.size() - 9];
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++; // 对snp开始计数
                    if (cnt % 1000000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }

                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();
                    List<String> l = new ArrayList<>();
                    l = PStringUtils.fastSplit(temp, "\t");
                    String chr = l.get(0);
                    String pos = l.get(1);
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
                    sb.append(String.format("%.4f", relativeMean));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println();
                System.out.println(f.getName() + "\t" + cnt + " sites is completed");
            } catch (Exception e) {
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }

        });
    }

    /**
     * chr001_5000.vcf --> chr001_5000_IndiHeter.txt Calculate the heterozygote
     * count and propotion by individual taxa
     */
    public void calIndiHeter(String infileDirS, String outfileDirS) {
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
                String infileS = f.getAbsolutePath();
                String indioutfileS = new File(outfileDirS, f.getName().split(".vc")[0] + "_IndiHeter.txt").getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                BufferedWriter indibw = IOUtils.getTextWriter(indioutfileS);
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
                 */
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
                indibw.flush();
                indibw.close();
                System.out.println(f.getName() + " is completed at " + indioutfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * chr001_5000.vcf --> chr001_5000_SNPheter.txt Calculate the heterozygote
     * count and propotion by individual taxa
     */
    public void calSNPSitesHeter(String infileDirS, String outfileDirS) {
        //String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/000_sampleVCF";
        //String outfileDirS = "/Users/Aoyue/Documents/test/";
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
            //Start to calwindow the time beginning
            //long startTime = System.nanoTime();
            //System.out.println("******************************************************");
            //Start to calwindow heterozygous
            try {
                String infileS = f.getAbsolutePath();
                String siteoutfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    siteoutfileS = new File(outfileDirS, f.getName().split(".vc")[0] + "_SNPheter.txt").getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    siteoutfileS = new File(outfileDirS, f.getName().split(".vc")[0] + "_SNPheter.txt").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(siteoutfileS);
                bw.write("Chr\tPos\tHetNum\tHomNum\tHetPropotion\n");
                //bw.write("Chr\tPos\tHetPropotion\n");
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
                 */
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
                        //
                        bw.write(te[0] + "\t" + te[1] + "\t" + String.format("%.0f", hetNum) + "\t" + String.format("%.0f", homNum) + "\t" + String.format("%.5f", hetRate) + "\n");
                        //bw.write(te[0]+ "\t"+ te[1] +  "\t"+ String.format("%.5f", hetRate) + "\n");
                    }
                }

                br.close();
                bw.flush();
                bw.close();
                System.out.println(f.getName() + " is completed at " + siteoutfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            //文件处理完毕，计时
//            long endTime = System.nanoTime();
//            float excTime = (float) (endTime - startTime) / 1000000000;
            //System.out.println("******************************************************" );
            //System.out.println("Execution time: " + String.format("%.2f", excTime) + "s" + "    or " + String.format("%.2f", excTime / 60) + " min");
        });
    }

    /**
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void calSNPHetMissMaf(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
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
                String outfileS = new File(outfileDirS, f.getName().split(".vc")[0] + "_SNPheterMissMaf.txt.gz").getAbsolutePath(); //输出是压缩的
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                //bw.write("Chr\tPos\tHetPropotion\n");
                //bw.write("Chr\tPos\tHetNum\tHetPropotion\tMissingNum\tMissProportion\tMaf\n");
                bw.write("Chr\tPos\tHetProportion\tMissProportion\tMaf");
                bw.newLine();
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
//                        bw.write(te[0] + "\t" + te[1] + "\t" + String.format("%.0f", hetNum) + "\t" + String.format("%.5f", hetRate)
//                                + "\t" + String.format("%.0f", missNum) + "\t" + String.format("%.5f", missRate) + "\t" + String.format("%.5f", maf) + "\n");
                        bw.write(te[0] + "\t" + te[1] + "\t" + String.format("%.5f", hetRate)
                                + "\t" + String.format("%.5f", missRate) + "\t" + String.format("%.5f", maf));
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
     * chr001.ABDgenome.filterMiss_subset.vcf.gz
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergesubsetVCF(String infileDirS, String outfileS) {
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
     * chr002.ABgenome.filterMiss.vcf --> chr002.ABgenome.filterMiss_subset.vcf
     * 过滤了alt含有2个allele的sites
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void subsetVCF(String infileDirS, String outfileDirS, String extractRatio) {



        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0] + "_bi.subset.vcf").getAbsolutePath(); //输出非压缩格式，到最后统一压缩
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnttotal++;
                        l = PStringUtils.fastSplit(temp);
                        if (l.get(4).contains(",") || (l.get(4).contains("D")) || (l.get(4).contains("I"))) {
                            continue; // 第4列是alt的信息，若有2个等位基因，则去除这一行
                        }
                        double r = Math.random();
                        double ratio = Double.parseDouble(extractRatio);

                        if (r > ratio) {
                            continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                        }

                        bw.write(temp);
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\twith " + cnttotal + " bp has a subset of\t" + cntsubset + "\tbiallelic SNPs is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     *
     * @param infileS
     * @param outfileS
     */
    public void subsetVCF_singleStream(String infileS, String outfileS, String ratio) {
        //String infileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001.vcf";
        //String outfileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001_subset.vcf";

        //String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/002_merge/subsetchr1_15.vcf.gz";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/004_pca/002_merge/subset10ksnp.vcf.gz";
        try {
            //BufferedReader br = IOUtils.getTextReader(infileS);
            //BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            String temp = null;
            int cnt = 0;
            System.out.println(new SimpleDateFormat().format(new Date()) + "    program execution.\n");
            long startTime = System.nanoTime();
            Double Ratio = Double.parseDouble(ratio);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                } else {
                    cnt++;
                    double r = Math.random();
                    if (r > Ratio) {
                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(4).contains(",")) {
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
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * chr005.Dlineage.vcf --> chr005.Dlineage.maf0.005.bi.vcf 过滤D和I，
     * 保留maf大于0.05的有1个或者2个alt的位点
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void filterIndelMaf(String infileDirS, String outfileDirS) {
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
        //System.out.println("FileName\tTotalSNPNum\tBiallelicMafmore0.01Num\tTriallelicMafmore0.01Num\tTriallelicAlt2more0.01Num\tProportionofTriallelicMafmore0.01Num\tProportionofTriallelicAlt2more0.01Num");
        System.out.println("Chr\tTotalSNP Num\tBiallelic Num(Maf>0.01)\tTriallelic Num(Maf>0.01)\tTriallelic Num(Alt2F>0.01)\tTriallelic Ratio(Maf>0.01)\tTriallelic Ratio(Alt2F>0.01)");
        fsList.parallelStream().forEach(f -> {
//        fsList.stream().forEach(f -> {  
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".maf0.01.SNP.vcf")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".maf0.01.SNP.vcf")).getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = null;
                String te[] = null;
                int cnt = 0;
                int biallelicMafmoreNum = 0;
                int cntCmaf12 = 0; //alt1 +alt2 大于0.005的个数
                int cntCmaf2 = 0; //alt2 大于 0.005的个数
                double ProportionofTriallelicMafmoreNum = Double.MIN_VALUE;
                double ProportionofTriallelicAlt2moreNum = Double.MIN_VALUE;

                while ((temp = br.readLine()) != null) {
                    int genoNum = 0;
                    double homNum = 0;
                    double hetNum = 0;
                    double hetRate = 0;
                    double missNum = 0;
                    double missRate = 0;

                    double refAlleleGametes = 0;
                    double alt1AlleleGametes = 0;
                    double alt2AlleleGametes = 0;
                    double refAF = 0;
                    double alt1AF = 0;
                    double alt2AF = 0;
                    double maf = 0;

                    if (temp.startsWith("#")) {//过滤含有#的注释部分
                        bw.write(temp);
                        bw.newLine();
                    } else { //开始进行位点的判断和计算
                        te = temp.split("\t");
                        String alt = PStringUtils.fastSplit(temp).get(4);
                        //***************************************************** 开始分情况讨论，这里分含逗号和不含逗号的情况 ***************************
                        if (te[4].length() == 1) { //不含有逗号的情况，即只有一个alt。又开始分，是D 是I 是ATGC 3种情况
                            if (alt.contains("D") || alt.contains("I")) {
                                continue; //只有一个alt且不是indel
                            }
                            cnt++;
                            for (int i = 9; i < te.length; i++) {
                                if (te[i].startsWith(".")) {
                                    missNum++;
                                }
                                if (!te[i].startsWith(".")) {
                                    genoNum++; //have the genotype
                                    if (te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                        hetNum++; //the number of heterozygous
                                        refAlleleGametes++;
                                        alt1AlleleGametes++;
                                    }
                                    if (te[i].startsWith("0/0")) {
                                        homNum++; //the number of heterozygous
                                        refAlleleGametes++;
                                        refAlleleGametes++;
                                    }
                                    if (te[i].startsWith("1/1")) {
                                        homNum++;
                                        alt1AlleleGametes++;
                                        alt1AlleleGametes++;
                                    }
                                }
                            }
                            //hetRate = hetNum / genoNum;
                            //missRate = missNum / (missNum + genoNum);
                            refAF = refAlleleGametes / (refAlleleGametes + alt1AlleleGametes);
                            alt1AF = alt1AlleleGametes / (refAlleleGametes + alt1AlleleGametes);;
                            if (refAF >= alt1AF) {
                                maf = alt1AF;
                            } else {
                                maf = refAF;
                            }

                            if (maf <= 0.01) {
                                continue;
                            }
                            biallelicMafmoreNum++;
                            StringBuilder sb = new StringBuilder();
                            bw.write(sb.append(temp).toString());
                            bw.newLine();

                        } else if (te[4].length() > 1) { //含有逗号的情况，即有2个alt。又开始分，是D 是I 是ATGC 3种情况
                            if (alt.contains("D") || alt.contains("I")) {
                                continue;
                            }
                            cnt++;
                            for (int i = 9; i < te.length; i++) {
                                if (te[i].startsWith(".")) {
                                    missNum++;
                                }
                                if (!te[i].startsWith(".")) {
                                    genoNum++; //have the genotype
                                    if (te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                        hetNum++; //the number of heterozygous
                                        refAlleleGametes++;
                                        alt1AlleleGametes++;
                                    }
                                    if (te[i].startsWith("0/0")) {
                                        homNum++; //the number of heterozygous
                                        refAlleleGametes++;
                                        refAlleleGametes++;
                                    }
                                    if (te[i].startsWith("1/1")) {
                                        homNum++;
                                        alt1AlleleGametes++;
                                        alt1AlleleGametes++;
                                    }
                                    if (te[i].startsWith("2/2")) {
                                        homNum++;
                                        alt2AlleleGametes++;
                                        alt2AlleleGametes++;
                                    }
                                    if (te[i].startsWith("0/2") || te[i].startsWith("2/0")) {
                                        hetNum++;
                                        refAlleleGametes++;
                                        alt2AlleleGametes++;
                                    }
                                    if (te[i].startsWith("1/2") || te[i].startsWith("2/1")) {
                                        hetNum++;
                                        alt1AlleleGametes++;
                                        alt2AlleleGametes++;
                                    }
                                }
                            }

                            //hetRate = hetNum / genoNum;
                            //missRate = missNum / (missNum + genoNum);
                            refAF = refAlleleGametes / (refAlleleGametes + alt1AlleleGametes + alt2AlleleGametes);
                            alt1AF = alt1AlleleGametes / (refAlleleGametes + alt1AlleleGametes + alt2AlleleGametes);
                            alt2AF = alt2AlleleGametes / (refAlleleGametes + alt1AlleleGametes + alt2AlleleGametes);

                            if (refAF >= alt1AF) { //先计算maf值大小，进行0.05的判断
                                maf = alt1AF;
                            } else {
                                maf = refAF;
                            }

                            if (maf <= 0.01) { //如果maf值小于0.005则过滤掉
                                continue;
                            }
                            cntCmaf12++;
                            if (alt2AF > 0.01) {
                                cntCmaf2++;
                            }
                            StringBuilder sb = new StringBuilder();
                            bw.write(sb.append(temp).toString());
                            bw.newLine();
                        }

                        //bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.5f", hetRate) + "\n");
                        //bw.write(te[0] + "\t" + te[1] + "\t" + String.format("%.0f", hetNum) + "\t" + String.format("%.5f", hetRate)
                        // + "\t" + String.format("%.0f", missNum) + "\t" + String.format("%.5f", missRate) + "\t" + String.format("%.5f", maf) + "\n");
                    }//////////////////////////在这里活动

                }
                ProportionofTriallelicMafmoreNum = (double) cntCmaf12 / (double) (biallelicMafmoreNum + cntCmaf12);
                ProportionofTriallelicAlt2moreNum = (double) cntCmaf2 / (double) (biallelicMafmoreNum + cntCmaf12);
                br.close();
                bw.flush();
                bw.close();

                System.out.println(String.valueOf(f.getName()) + "\t" + cnt + "\t" + String.valueOf(biallelicMafmoreNum) + "\t" + cntCmaf12 + "\t" + cntCmaf2 + "\t" + ProportionofTriallelicMafmoreNum + "\t" + ProportionofTriallelicAlt2moreNum + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * chr005.Dlineage.vcf --> chr005.Dlineage.maf0.005.bi.vcf 过滤D I
     * 和含有3个等位位点的pos，保留只有一个alt并且maf大于0.005的pos!!!保留只有一个alt并且maf大于0.005的pos.
     *
     * @param infileDirS
     */
    public void filterAlleleMaf(String infileDirS, String outfileDirS) {
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
        System.out.println("FileName\tbiallelicNum\tBiallelicMafmore0.005Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".maf0.005.bi.vcf")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".maf0.005.bi.vcf")).getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = null;
                String te[] = null;
                int biallelicNum = 0;
                int biallelicMafmoreNum = 0;
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
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        te = temp.split("\t");
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
                            //hetRate = hetNum / genoNum;
                            //missRate = missNum / (missNum + genoNum);
                            refAF = refAlleleGametes / (refAlleleGametes + altAlleleGametes);
                            altAF = altAlleleGametes / (refAlleleGametes + altAlleleGametes);;
                            if (refAF >= altAF) {
                                maf = altAF;
                            } else {
                                maf = refAF;
                            }

                            if (maf <= 0.01) {
                                continue;
                            }
                            biallelicMafmoreNum++;
                            StringBuilder sb = new StringBuilder();
                            bw.write(sb.append(temp).toString());
                            bw.newLine();

                        }
                    }
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

    /**
     * chr005.Dgenome.vcf.gz --> chr005.Dgenome.bi.vcf.gz Keep only the binary
     * allele 过滤VCF文件，只保留双等位位点，即含有1个alt
     *
     * @param infileS
     * @param outfileS
     */
    public void filterSNPtoBi_singleThread(String infileS, String outfileS) {
        File f = new File(infileS);
        try {
            String chr = f.getName().substring(3, 6); //提取染色体号 001
            int chrint = Integer.parseInt(chr); //将染色体号转化为数字
            BufferedReader br = AoFile.readFile(f.getAbsolutePath());
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            int biallelicNum = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                } else {
                    cnt++;
                    temp = temp.substring(0,50);
                    String alt = PStringUtils.fastSplit(temp).get(4);
                    if (alt.contains(",")) continue;
                    if (alt.equals("D")) continue;
                    if (alt.equals("I")) continue;
                    biallelicNum++;
                    bw.write(temp);
                    bw.newLine();
                }

            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(chrint + "\t" + biallelicNum + "\tBi-SNP\t" + cnt + "\tTotalVariants\t" + f.getName() + " is completed at " + outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * chr005.Dgenome.vcf.gz --> chr005.Dgenome.bi.vcf.gz Keep only the binary
     * allele 过滤VCF文件，只保留双等位位点，即含有1个alt
     *
     * @param infileDirS
     */
    public void filterSNPtoBi_parallel(String infileDirS, String outfileDirS) {
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        File[] fs1 = IOUtils.listFilesEndsWith(fs,"vcf.gz");
        File[] fs2 = IOUtils.listFilesEndsWith(fs,"vcf");
        List<File> fsList = new ArrayList<>();
        for (int i = 0; i < fs1.length; i++) {
            fsList.add(fs1[i]);
        }
        for (int i = 0; i < fs2.length; i++) {
            fsList.add(fs2[i]);
        }

        Collections.sort(fsList);
        //        System.out.println("Chr\tSNPNum\tBiallelicNum\tIndelNum\tInsertionNum\tDelectionNum\t");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                String outfileS = new File(outfileDirS, "chr" + chr + "_Bi.vcf").getAbsolutePath();
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
//                int snpNum = 0;
                int biallelicNum = 0;
//                int indelNum = 0;
//                int insertionNum = 0;
//                int delectionNum = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {

                        String alt = PStringUtils.fastSplit(temp).get(4);
//                        if (alt.contains("D")) {
//                            delectionNum++;
//                        }
//                        if (alt.contains("I")) {
//                            insertionNum++;
//                        }

//                        if (!(alt.contains(",")) && !(alt.equals("D")) && !(alt.equals("I"))) {
//                            biallelicNum++;
//                            bw.write(temp);
//                            bw.newLine();
//                        }

                        if (alt.contains(",")) continue;
                        if (alt.equals("D")) continue;
                        if (alt.equals("I")) continue;
                        biallelicNum++;
                        bw.write(temp);
                        bw.newLine();
                    }

                }
//                indelNum = delectionNum + insertionNum;
//                snpNum = cnt - indelNum;
                br.close();
                bw.flush();
                bw.close();
//                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(snpNum) + "\t" + String.valueOf(biallelicNum) + "\t" + String.valueOf(indelNum) + "\t" + String.valueOf(insertionNum) + "\t" + String.valueOf(delectionNum));
                System.out.println(chrint + "\t" + biallelicNum + "\t" + f.getName() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 统计一下位点既含有D又含有I的情况,确实存在着这种情况！！！
     * 此方法没有考虑既含有D又含有I的情况
     *
     * @param infileDirS
     * @deprecated
     */
    public void countrepeatIndelinFastCallformat(String infileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);

        System.out.println("Chr\tSNPNum\tRepeatDI\tBiallelicNum\tIndelNum\tInsertionNum\tDelectionNum\t");
        fsList.stream().forEach(f -> {
            try {
                if (f.getName().contains("036")) {
                    String infileS = f.getAbsolutePath();
                    BufferedReader br = null;
                    if (infileS.endsWith(".vcf")) {
                        br = IOUtils.getTextReader(infileS);
                    } else if (infileS.endsWith(".vcf.gz")) {
                        br = IOUtils.getTextGzipReader(infileS);
                    }
                    String chr = f.getName().substring(3, 6); //提取染色体号 001
                    int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                    String temp = null;
                    int cnt = 0;
                    int snpNum = 0;
                    int biallelicNum = 0;
                    int indelNum = 0;
                    int insertionNum = 0;
                    int delectionNum = 0;
                    int repeatDI = 0;
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
                        if (alt.contains("I") && alt.contains("D")) {
                            repeatDI++;
                        }

                        if (!(alt.contains(",")) && !(alt.equals("D")) && !(alt.equals("I"))) {
                            biallelicNum++;
                        }
                    }

                    indelNum = delectionNum + insertionNum;
                    snpNum = cnt - indelNum; //只要含有D 和 I， 就不算是SNP
                    br.close();
                    System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(snpNum) + "\t" + repeatDI + "\t" + String.valueOf(biallelicNum) + "\t" + String.valueOf(indelNum) + "\t" + String.valueOf(insertionNum) + "\t" + String.valueOf(delectionNum));

                }

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    /**
     * 对已生成的14条染色体进行计数.注意加log文件，结果在log文件中显示。
     *  最终采用方法。
     *
     * @param infileDirS
     */
    public static void countSitesinFastCallformat_fromTxt(String infileDirS) {

        List<File> fsList = AoFile.getFileListInDir(infileDirS);

        System.out.println("Chr\tN_RawSNPs\tN_BiallelicSNPs\tN_TriallelicSNPs\tIndels\tInsertions\tDeletions");
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                int cntSNP = 0;
                int cntBi = 0;
                int cntTri = 0;
                int cntIndel = 0;
                int cntI = 0;
                int cntD = 0;
                String temp = null;
                String header = br.readLine();
                while ((temp = br.readLine()) != null) { //是否含有D I
                    String alt = PStringUtils.fastSplit(temp).get(3);
                    /**
                     * 含有逗号，即有2个alt.
                     * case 1: A,T
                     * case 2: D,I
                     * case 3: D,G
                     * case 4: I,T
                     */
                    if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                        boolean ifD = false;
                        if (!alt.contains("D") && (!alt.contains("I"))) {
                            cntTri++;
                            cntSNP++;
                        }
                        if (alt.contains("D")) {
                            cntD++;
                            cntIndel++;
                            ifD = true;
                        }
                        if (alt.contains("I")) {
                            cntI++;
                            if (ifD == false) { //针对 case 2的情况，即含有D又含有I， 这时在上面的D中已经加过 cntIndel了，所以不用再加了
                                cntIndel++;
                            }
                        }

                    } else if (alt.length() == 1) { //1个alt的情况;
                        if (!alt.equals("D") && (!alt.equals("I"))) {
                            cntBi++;
                            cntSNP++;
                        }
                        if (alt.equals("D")) {
                            cntD++;
                            cntIndel++;
                        }
                        if (alt.equals("I")) {
                            cntI++;
                            cntIndel++;
                        }
                    }
                }
                br.close();
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cntSNP) + "\t" + String.valueOf(cntBi) + "\t" + String.valueOf(cntTri) + "\t" + String.valueOf(cntIndel) + "\t" + String.valueOf(cntI) + "\t" + String.valueOf(cntD));
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
    public static void countSitesinFastCallformat(String infileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        File[] fs1 = IOUtils.listFilesEndsWith(fs,"vcf.gz");
        File[] fs2 = IOUtils.listFilesEndsWith(fs,"vcf");
        List<File> fsList = new ArrayList<>();
        for (int i = 0; i < fs1.length; i++) {
            fsList.add(fs1[i]);
        }
        for (int i = 0; i < fs2.length; i++) {
            fsList.add(fs2[i]);
        }

        Collections.sort(fsList);
        System.out.println("Chr\tN_RawSNPs\tN_BiallelicSNPs\tN_TriallelicSNPs\tIndels\tInsertions\tDeletions");
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                int cntSNP = 0;
                int cntBi = 0;
                int cntTri = 0;
                int cntIndel = 0;
                int cntI = 0;
                int cntD = 0;
                String temp = null;
                while ((temp = br.readLine()) != null) { //是否含有D I
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temp = temp.substring(0, 150);
                    String alt = PStringUtils.fastSplit(temp).get(4);
                    /**
                     * 含有逗号，即有2个alt.
                     * case 1: A,T
                     * case 2: D,I
                     * case 3: D,G
                     * case 4: I,T
                     */
                    if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                        boolean ifD = false;
                        if (!alt.contains("D") && (!alt.contains("I"))) {
                            cntTri++;
                            cntSNP++;
                        }
                        if (alt.contains("D")) {
                            cntD++;
                            cntIndel++;
                            ifD = true;
                        }
                        if (alt.contains("I")) {
                            cntI++;
                            if (ifD == false) { //针对 case 2的情况，即含有D又含有I， 这时在上面的D中已经加过 cntIndel了，所以不用再加了
                                cntIndel++;
                            }
                        }

                    } else if (alt.length() == 1) { //1个alt的情况;
                        if (!alt.equals("D") && (!alt.equals("I"))) {
                            cntBi++;
                            cntSNP++;
                        }
                        if (alt.equals("D")) {
                            cntD++;
                            cntIndel++;
                        }
                        if (alt.equals("I")) {
                            cntI++;
                            cntIndel++;
                        }
                    }
                }
                br.close();
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cntSNP) + "\t" + String.valueOf(cntBi) + "\t" + String.valueOf(cntTri) + "\t" + String.valueOf(cntIndel) + "\t" + String.valueOf(cntI) + "\t" + String.valueOf(cntD));
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 根据1 -42 条染色体的统计，转换成 1A- 7A 形式的统计，再装换成 A sub  B sub 的统计；
     * 本方法只解决 1A 2A 3A 到 A sub  的统计
     */
    public static void mergeChr1Aand2A_bysubgenome(String infileS, String outfileS){
        int cntcolumn = AoFile.countFileColumnNumber(infileS)-1;
        File parent = new File(outfileS).getParentFile();
        String outfileS1 = new File(parent,"Asub_aoyue.txt").getAbsolutePath();
        String outfileS2 = new File(parent,"Bsub_aoyue.txt").getAbsolutePath();
        String outfileS3 = new File(parent,"Dsub_aoyue.txt").getAbsolutePath();

        try {
            List<Integer>[] valueArray = new List[cntcolumn];
            for (int i = 0; i < valueArray.length; i++) {
                valueArray[i] = new ArrayList<>();
            }
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS1);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0).substring(1);
                if (!chr.equals("A")) continue;
                for (int i = 1; i < l.size(); i++) {
                    int index = i - 1;
                    valueArray[index].add(Integer.parseInt(l.get(i)));
                }
                bw.write(temp);
                bw.newLine();
            }

            bw.write("Asub");
            for (int i = 0; i < valueArray.length; i++) {
                int sum = AoMath.listSum_byint(valueArray[i]);
                bw.write("\t" + sum);
            }
            bw.newLine();
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            List<Integer>[] valueArray = new List[cntcolumn];
            for (int i = 0; i < valueArray.length; i++) {
                valueArray[i] = new ArrayList<>();
            }
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS2);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0).substring(1);
                if (!chr.equals("B")) continue;
                for (int i = 1; i < l.size(); i++) {
                    int index = i - 1;
                    valueArray[index].add(Integer.parseInt(l.get(i)));
                }
                bw.write(temp);
                bw.newLine();
            }

            bw.write("Bsub");
            for (int i = 0; i < valueArray.length; i++) {
                int sum = AoMath.listSum_byint(valueArray[i]);
                bw.write("\t" + sum);
            }
            bw.newLine();
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            List<Integer>[] valueArray = new List[cntcolumn];
            for (int i = 0; i < valueArray.length; i++) {
                valueArray[i] = new ArrayList<>();
            }
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS3);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0).substring(1);
                if (!chr.equals("D")) continue;
                for (int i = 1; i < l.size(); i++) {
                    int index = i - 1;
                    valueArray[index].add(Integer.parseInt(l.get(i)));
                }
                bw.write(temp);
                bw.newLine();
            }

            bw.write("Dsub");
            for (int i = 0; i < valueArray.length; i++) {
                int sum = AoMath.listSum_byint(valueArray[i]);
                bw.write("\t" + sum);
            }
            bw.newLine();
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        File f1 = new File(outfileS1);
        File f2 = new File(outfileS2);
        File f3 = new File(outfileS3);

        File[] fs = {f1,f2,f3};
        AoFile.mergeTxt_byFileArray(fs,outfileS);
    }

    /**
     * 将计算出的snp位点数进行合并
     */
    public void mergeChr1and2_ABgenome(String infileS, String outfileS) {
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/002_countSites/countSites_mergeChr1and2.txt";

        //infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK.txt";
        //outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/004_gVCF/002_vcf_GATK/002_countSites/countSites_fromGATK_mergeChr1and2.txt";
        String[] chr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15,16,17,18,19,20,21};
        HashMap<Integer, String> hmcntchr = new HashMap<>();
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
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int[] site1 = new int[l.size()];
                int[] site2 = new int[l.size()];
                int[] site = new int[l.size()];
                HashMap<Integer, Integer>[] hmcnt = new HashMap[l.size()];
                for (int i = 1; i < l.size(); i++) {
                    site1[i] = Integer.parseInt(l.get(i));
                }

                if ((temp = br.readLine()) != null) { //then the second line
                    cnt++;
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 1; i < l.size(); i++) { //record site2
                        site2[i] = Integer.parseInt(l.get(i));
                    }
                    for (int i = 1; i < l.size(); i++) { //calculate the sum site
                        site[i] = site1[i] + site2[i];
                    }

                    bw.write(hmcntchr.get(cnt)); // + "\t" + hmcnt.get(cnt));
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + site[i]);
                    }
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
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int[] site1 = new int[l.size()];
                int[] site2 = new int[l.size()];
                int[] site = new int[l.size()];
                HashMap<Integer, Integer>[] hmcnt = new HashMap[l.size()];
                for (int i = 1; i < l.size(); i++) {
                    site1[i] = Integer.parseInt(l.get(i));
                }

                if ((temp = br.readLine()) != null) { //then the second line
                    cnt++;
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 1; i < l.size(); i++) { //record site2
                        site2[i] = Integer.parseInt(l.get(i));
                    }
                    for (int i = 1; i < l.size(); i++) { //calculate the sum site
                        site[i] = site1[i] + site2[i];
                    }

                    bw.write(hmcntchr.get(cnt)); // + "\t" + hmcnt.get(cnt));
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + site[i]);
                    }
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
     * 针对特定的小麦42条染色体，
     * 将表格汇总的每一列数值按照1A 1B 1D 相加，产生一个新的表格，表格中可以有很多列的值,但必须为int类型
     *
     * @param infileS
     * @param outfileS
     */
    public static void  mergeChr1and2txt_int(String infileS, String outfileS) {
        String[] chr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
        HashMap<Integer, String> hmcntchr = new HashMap<>();

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

            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int[] site1 = new int[l.size()];
                int[] site2 = new int[l.size()];
                int[] site = new int[l.size()];
                HashMap<Integer, Integer>[] hmcnt = new HashMap[l.size()];
                for (int i = 1; i < l.size(); i++) {
                    site1[i] = Integer.parseInt(l.get(i));
                }

                if ((temp = br.readLine()) != null) { //then the second line
                    cnt++;
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 1; i < l.size(); i++) { //record site2
                        site2[i] = Integer.parseInt(l.get(i));
                    }
                    for (int i = 1; i < l.size(); i++) { //calculate the sum site
                        site[i] = site1[i] + site2[i];
                    }

                    bw.write(hmcntchr.get(cnt)); // + "\t" + hmcnt.get(cnt));
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + site[i]);
                    }
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
     * 将表格汇总的每一列数值按照1A 1B 1D 相加，产生一个新的表格，表格中可以有很多列的值,为int类型 或者 double 类型
     *
     * @param infileS
     * @param outfileS
     */
    public void mergeChr1and2txt_double(String infileS, String outfileS) {
        String[] chr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int[] cnts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
        HashMap<Integer, String> hmcntchr = new HashMap<>();

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

            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                Double[] site1 = new Double[l.size()];
                Double[] site2 = new Double[l.size()];
                Double[] site = new Double[l.size()];
                HashMap<Integer, Integer>[] hmcnt = new HashMap[l.size()];
                for (int i = 1; i < l.size(); i++) {
                    site1[i] = Double.parseDouble(l.get(i));
                }

                if ((temp = br.readLine()) != null) {
                    cnt++;
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 1; i < l.size(); i++) {
                        site2[i] = Double.parseDouble(l.get(i));
                    }
                    for (int i = 1; i < l.size(); i++) {
                        site[i] = site1[i] + site2[i];
                    }

                    bw.write(hmcntchr.get(cnt)); // + "\t" + hmcnt.get(cnt));
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + site[i]);
                    }
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
     * 将计算出的snp位点数进行合并，成1A 1B 1D形式；
     *
     */
    @Deprecated
    public void mergeChr1and2_deprecated(String infileS, String outfileS) {
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
     */
    public static int countSites_fromVCF(String infileS) {
        int cnt = 0;
        try {
            BufferedReader br = AoFile.readFile(infileS);
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
        System.out.println("Total SNP num is "+ String.valueOf(cnt));

        return cnt;
    }

    /**
     * Count the snp sites via single stream and print it into the inputstream.
     *
     * @param infileDirS
     */
    public static void countSites_singleStream(String infileDirS) {
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = 0;
            try {
                String infileS = fs[i].getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
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
        System.out.println("Total SNP num is "+ String.valueOf(sum));
    }
}
