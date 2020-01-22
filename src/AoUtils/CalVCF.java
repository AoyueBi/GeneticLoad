/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import format.table.RowTable;
import utils.IOUtils;
import utils.PStringUtils;

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
            BufferedReader br = null;
            BufferedWriter bw = null;
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
            System.out.println(infileS + " is completed at " + outfileS + "\tActual taxa size: " + indexHexa.size() + "\tGoal taxa size : " + lhexa.size());
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
    public boolean ifSegregation(String[] genoArray, String altList) {
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
     * return the site heterozygosity from vcf
     *
     * @param genoArray
     * @return
     */
    public Double calSNPSitesHeter(String[] genoArray){
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
