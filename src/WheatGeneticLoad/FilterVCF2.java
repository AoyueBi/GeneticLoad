package WheatGeneticLoad;

import AoUtils.*;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;

import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;

public class FilterVCF2 {

    public FilterVCF2(){

        /**
         * 过滤 MAF > 0.01  Occurrence>2 MissingRate<0.2
         */
        this.filter_parallel(); //老师的方法
//        this.filterMafbyPopHTD(); //我的方法
//        String a = "";
//        String b = "";
//        this.filter_singleThread(a,b); //最终过滤的时候采用的方法
//        this.script();
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_script/sh_filterMAFmissOccurrence20200522.sh",3,11);
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_script/sh_filterMAFmissOccurrence_2_20200522.sh",4,2);
        // 统计 Indel 位点个数
//        CountSites.mergeChr1and2txt_int("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/log_043_countSitesinFastCallformat_fixVMap2.0Indels_20200818.txt","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/CountVariants_fixVMap2.0_Indel_20200818.txt");


        /**
         * 对过滤MAF miss occurrence 的程序进行验证
          */
//        this.extractField();
//        this.getOccurrenceVCF();
//        this.getOccurrenceByTaxa();
//        this.identifyFinalResult();

        /**
         * fix VMap2.0
         */
//        this.modifyVMap2(); //修改名字
//        new CountSites().mergeChr1and2txt_int("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/log_043_countSitesinFastCallformat_fixVMap2.0_20200522.txt","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/CountVariants_fixVMap2.0_20200522.txt");
//        CountSites.mergeChr1and2txt_int("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/log_043_countSitesinFastCallformat_fixVMap2.0_20200601.txt","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/CountVariants_fixVMap2.0_202006.txt");
//        CountSites.mergeChr1and2txt_int("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/log_043_countSitesinFastCallformat_fixVMap2.0_20200604.txt","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/CountVariants_fixVMap2.0_202006.txt");

//        this.bgzip();
//        this.sortTaxaName();

        /**
         *  VCF quality control
         */

//        this.mergeVCF();
//        this.QC();
//        this.mergeCheckFile();
//        this.getBinTable();
//        this.statVcfDepth_SD(); //计算按照亚群和按照倍性区分的深度
//        this.mergeTxtandAddGroup();

//        this.mkDepthOfVMapII(); //计算taxa的深度
//        this.mkDepthSummary();
//        this.mergeTaxaDepth();
//        this.calSite();

//        this.getMergedSubsetVCF_Hexaploid();

        /**
         * make vmap2.1
         */

//        this.scriptFilterSNPtoBi();

        /**
         * 提取六四二倍体的VCF
         */

//        this.runJarParallele();
//        SplitScript.splitScript2("/Users/Aoyue/Documents/sh_vmap2.0tovmap2.1_20200526.sh",21,2);
//        this.getsharedSNP(); //获取六倍体四倍体二倍体共有的SNP
//        this.mergeSharedSNP();
//        CountSites.mergeChr1Aand2A_bysubgenome("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/003_pop/001_getsharedSNP_merged.txt","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/003_pop/001_getsharedSNP_merged_bysub.txt");


        /**
         * Variants rate
         */

//        this.mergeVariantRate();

        /**
         * ************** Indel的质控, correct 前后变化
         */
//        this.getIndelVCF();
//        this.extractPosAllele();

//        this.QC_indel();
//        this.mergeCheckFile();
//        this.getBinTable();
//        this.statVcfDepth_SD();
//        this.mergeTxtandAddGroup();
        /**
         * ************** VMap2.0 finalize
         */

//        this.finalizeVMapII();
//        CountSites.mergeChr1Aand2A_bysubgenome("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/CountVariants_fixVMap2.0_Indel_20200818.txt","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/002_vmap2.0/VMap2.0_Inde_bySubgenome.txt");

//        this.test();

    }


    public void test(){

        for (int i = 1; i < 43; i++) {
            if (i != 29)continue;
            System.out.println(i);
        }
    }

    /**
     * 将2020-05 VMap2.0 (i.e. hapscanner并过滤 MAF>=0.01 Miss<=0.2 Occurrence >=2) 中的 SNP 和Indel 结果只保留SNP(include biallelic and tri-allelic SNP)
     * 将2020-08 新做出来Indel (revised hapscanner output --> 并过滤 MAF>=0.01 Miss<=0.2 Occurrence >=2) 和上文SNP排序并合并
     */
    public void finalizeVMapII(){

//        String infileDirS = ""; //VMap2.0 include SNP and error Indel
//        String infile2DirS = ""; //Vmap2.0 only include Indel
//        String outfileDirS = "";


        //test data
//        String infileDirS = "/Users/Aoyue/Documents/in1"; //VMap2.0 include SNP and error Indel
//        String infile2DirS = "/Users/Aoyue/Documents/in2"; //Vmap2.0 only include Indel
//        String outfileDirS = "/Users/Aoyue/Documents/out";

        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_VMap2.0"; //VMap2.0 include SNP and error Indel
        String infile2DirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/107_VMap2.0_Indel"; //Vmap2.0 only include Indel
        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/108_VMap2.0_final";

//        int[] chrArray = new int[27];
//        for (int i = 2; i < 29; i++) {
//            chrArray[i-2]=i;
//        }
//        Arrays.sort(chrArray);

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        for (int i = 0; i < fs.length; i++) {
            int chr = Integer.parseInt(fs[i].getName().substring(3,6));
//            int id = Arrays.binarySearch(chrArray,chr);
//            if (id>-1)continue;

            if (chr != 29) continue;
            String infileS = fs[i].getAbsolutePath();
            String infileS2 = new File(infile2DirS,fs[i].getName().replaceFirst(".vcf",".vcf.gz")).getAbsolutePath();
            String outfileS = new File(outfileDirS,fs[i].getName()).getAbsolutePath();
            List<Record> rl = new ArrayList<>();
            try{
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);

                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt1 = 0;
                while ((temp = br.readLine()) != null){
                    if (temp.startsWith("#")){
                        bw.write(temp);bw.newLine();continue;
                    }
                    l = PStringUtils.fastSplit(temp);
                    String alt = l.get(4);
                    if (alt.contains("D")|alt.contains("I"))continue;
                    Record r = new Record(Integer.parseInt(l.get(1)),temp);
                    rl.add(r);
                    cnt1++;
                }
                br.close();

                System.out.println("chr " + chr + "\t******** record list 1 is completed with snp number\t" + cnt1);
                System.out.println(new SimpleDateFormat().format(new Date()));


                br = AoFile.readFile(infileS2);
                int cnt2 = 0;
                while ((temp = br.readLine()) != null){
                    if (temp.startsWith("#"))continue;
                    l = PStringUtils.fastSplit(temp);
                    Record r = new Record(Integer.parseInt(l.get(1)),temp);
                    rl.add(r);
                    cnt2++;
                }
                br.close();
                System.out.println("chr " + chr + "\t******** record list 2 is completed with indel number\t" + cnt2);
                System.out.println(new SimpleDateFormat().format(new Date()));

                Collections.sort(rl);
                System.out.println("chr " + chr + "\t******** finish sort the pos");
                System.out.println(new SimpleDateFormat().format(new Date()));
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < rl.size(); j++) {
                    sb.setLength(0);
                    sb.append(rl.get(j).r);
                    bw.write(sb.toString());bw.newLine();
                }
                bw.flush();bw.close();
                System.out.println("chr " + chr + "\t======== finish writing the file");
                System.out.println(new SimpleDateFormat().format(new Date()));


            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
        //nohup java -Xms400g -Xmx400g -jar PlantGenetics.jar > log_finalizeVMapII_20200820.txt 2>&1 &
        // cat /data4/home/aoyue/vmap2/aaPlantGenetics/log_finalizeVMapII_20200820.txt

        //nohup java -Xms400g -Xmx400g -jar PlantGenetics.jar > log_finalizeVMapII_20200820_chr29.txt 2>&1 &
        // cat /data4/home/aoyue/vmap2/aaPlantGenetics/log_finalizeVMapII_20200820_chr29.txt


    }

    class Record implements Comparable<Record>{
        int pos;
        public String r;
        public Record (int pos, String r){
            this.pos = pos;
            this.r = r;
        }
        @Override
        public int compareTo(Record o){return this.pos-o.pos;}
    }

    public void QC_indel(){
//        String abdinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/000_VCFwithonlyIndel/ABsubgenome_hexa.vcf.gz";
//        String abinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/000_VCFwithonlyIndel/ABsubgenome_tetra.vcf.gz";
//        String dinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/000_VCFwithonlyIndel/Dsubgenome_diploid.vcf.gz";
//        String abd_DsubinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/000_VCFwithonlyIndel/Dsubgenome_hexa.vcf.gz";

//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/001_ori";

//        this.checkQuality(abinfileS,outfileDirS,"AB");
//        this.checkQuality(abdinfileS,outfileDirS,"ABD");
//        this.checkQuality(dinfileS,outfileDirS,"D");
//        this.checkQuality(abd_DsubinfileS,outfileDirS,"ABD_Dsub");


/**
 * 抽取 chr036 比较修改程序前后的基因型变化，对maf造成的影响
 */

///////******** corrected Indel
//        String abdinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/chr036.vcf.gz";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/";
//        this.checkQuality(abdinfileS,outfileDirS,"ABD");

        //maf的分bin
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/ABD_site_QC.txt.gz";
//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.009;
//        double max = 0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);

///////******** incorrected Indel
//        String abdinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/errorFile/chr036_Indel.vcf.gz";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/errorFile/";
//        this.checkQuality(abdinfileS,outfileDirS,"ABD");


//        //maf的分bin
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/errorFile/ABD_site_QC.txt.gz";
//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.009;
//        double max = 0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/errorFile/maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);

        ///////******** corrected Insertion
//        String abd_Insertion_infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/Insertion.vcf.gz";
//        String abd_Deletion_infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/Deletion.vcf.gz";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/";
//        this.checkQuality(abd_Insertion_infileS,outfileDirS,"ABD_Insertion");
//        this.checkQuality(abd_Deletion_infileS,outfileDirS,"ABD_Deletion");

        //maf的分bin
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_Insertion_site_QC.txt.gz";
//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.009;
//        double max = 0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_insertion_maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);


//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_Deletion_site_QC.txt.gz";
//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.009;
//        double max = 0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_deletion_maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_insertion_maf.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_deletion_maf.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/036_hapScanner_Indel/test/insertion_deletion/ABD_insertion_deletion_maf.txt";
        File f = new File(infileS);
        File f2 = new File(infileS2);
        File[] fs = {f,f2};
        AoFile.mergeTxt_byFileArray(fs, outfileS);

    }

    /**
     * 提取Indel的chr pos ref alt
     */
    public void extractPosAllele(){
//        String infileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/001_IndelVCF";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/003_variationLibrary_Indel";
//        System.out.println("java -jar 007_mkHapPosAllele.jar " + infileDirS + " " + outfileDirS + " > log_mkHapPosAllele_20200816.txt 2>&1 &");

        String infileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/001_IndelVCF";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/004_hapPos";
        CalVCF.mkHapPosWithoutHeader(infileDirS,outfileDirS);

        //java -jar PlantGenetics.jar > log_getHapPos_20200816.txt 2 >&1 &

    }


    /**
     * 由于 hapscanner后，Indel的数目是 14,213,456， 故决定将indel全部抽取出来，进行质控
     * 质控内容包括 maf miss
     *
     */
    public void getIndelVCF(){
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/005_IndelsFrom_rawMergedVCF/001_IndelVCF";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0] + "_Indel.vcf").getAbsolutePath(); //输出非压缩格式，到最后统一压缩
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
                        if ((l.get(4).contains("D")) || (l.get(4).contains("I"))) {
                            bw.write(temp);
                            bw.newLine();
                            cntsubset++;
                        }

                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\twith " + cnttotal + " bp has\t" + cntsubset + "\tindels at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        // java -jar 050_getIndelVCF.jar > log_050_getIndelVCF_20200815.txt 2>&1 &
    }

    public void mergeVariantRate(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/017_variantsRate/001_variantsRate.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/017_variantsRate/002_variantsRate_byRefChr.txt";
//        new CountSites().mergeChr1and2txt_double(infileS,outfileS);

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/017_variantsRate/002_variantsRate_byRefChr.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/017_variantsRate/003_variantsRate_bySub.txt";
        CountSites.mergeChr1Aand2A_bysubgenome(infileS,outfileS);
    }


    /**
     * 从抽样的VCF中，提取六倍体的42条染色体的VCF文件，并合并成一个文件
     *
     */
    public void getMergedSubsetVCF_Hexaploid(){ //一次性将所有的jar都运行上


        //#/****************************** hexaploid ******************************/#
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/001_fromVMap2.0_singleChr0.001";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/025_subsetVCF/003_hexaploid_singleChr";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/log";
        String taxaListS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/BreadWheat_S420.txt";

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1.vcf").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_hexaploid.vcf").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            System.out.println("java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaListS + " > " + logfileS + " 2>&1 &");
        }

        new CountSites().mergesubsetVCF(outfileDirS,"/data4/home/aoyue/vmap2/analysis/025_subsetVCF/002_mergeVCFtoSub/ABDsubgenome_hexa.vcf.gz");
    }


    public void mergeSharedSNP(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/003_pop/log_getsharedSNP_20200621.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/029_countSiteSummary/003_pop/001_getsharedSNP_merged.txt";
        new CountSites().mergeChr1and2txt_int(infileS,outfileS);

    }

    /**
     * 输出ABD和AB共有的SNP位置
     */
    public void getsharedSNP(){
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/104_VCFbyPop/001_byPloidy/003_hexaploid";
        String abfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/104_VCFbyPop/001_byPloidy/002_tetraploid";
        String dfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/104_VCFbyPop/001_byPloidy/001_diploid";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/026_sharedSNPpos";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fsList);
        System.out.println("Chr" + "\t" + "Num_total_AABBDD" + "\t" + "Num_unique_AABBDD" + "\t" + "Num_shared" + "\t" + "Num_unique_AABBorDD" + "\t" + "Num_tatol_AABBorDD");
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String chr = f.getName().substring(3,6);
            int chrID = Integer.parseInt(chr);
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chrID);
            String infileDirS2 = null;
            String infileS2 = null;
            String outfileS = new File(outfileDirS,f.getName().substring(0,6) + "_sharedSNP_pos.txt").getAbsolutePath();
            if (subgenome.equals("D")) {
                infileDirS2 = dfileDirS;
                infileS2 = new File(infileDirS2,"chr" + chr + "_vmap2.1_diploid.vcf.gz").getAbsolutePath();

            }else{
                infileDirS2 = abfileDirS;
                infileS2 = new File(infileDirS2,"chr" + chr + "_vmap2.1_tetraploid.vcf.gz").getAbsolutePath();
            }

            int cntHexaSNP = 0;
            int cntsharedSNP = 0;
            int cntuniqueHexa = 0;
            int cntTDSNP = 0;
            int cntuniqueTD = 0;

            TIntArrayList polistHexa = CalVCF.extractVCFPos(infileS);
            polistHexa.sort();

            try {
                BufferedReader br = AoFile.readFile(infileS2);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write("Pos");
                bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    temp = temp.substring(0,20);
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    cntTDSNP++;
                    int index = polistHexa.binarySearch(pos);
                    if (index > -1){ //说明搜到了，是共享的
                        cntsharedSNP++;
                        bw.write(l.get(1));
                        bw.newLine();
                    }else if (index < 0){ //说明没搜到，是独有的
                        cntuniqueTD++;
                    }
                }

                cntHexaSNP = polistHexa.size();
                cntuniqueHexa = cntHexaSNP - cntsharedSNP;
                System.out.println(chrID + "\t" + cntHexaSNP + "\t" + cntuniqueHexa + "\t" + cntsharedSNP + "\t" + cntuniqueTD + "\t" + cntTDSNP);
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            // java -Xmx200g -jar PlantGenetics.jar > log_getsharedSNP_20200621.txt 2>&1 &
        });
    }







    //java -Xms50g -Xmx200g -jar 028_extractVCF.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII/chr001_vmap2.1.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/012_VCFbyPop/001_byPloid/hexaploid/chr001_vmap2.1_hexaploid.vcf /data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt > log_028/log_extractVCF_chr001_hexaploid20191107.txt 2>&1 &

    public void runJarParallele(){ //一次性将所有的jar都运行上

        // java -jar PlantGenetics.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/102_VMap2.0 /data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1 > log_filterAlleletoBi_20200525.txt 2>&1 &

        //#/****************************** hexaploid ******************************/#
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1";
//        String outfileDirS ="/data4/home/aoyue/vmap2/genotype/mergedVCF/104_VCFbyPop/001_byPloidy";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20200526";
//        String taxaListS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/BreadWheat_S420.txt";

//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/108_VMap2.0_final";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/032_vcf_S200_hexaploid";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20200921";
//        String taxaListS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/103_S196_hexaploid_10Xcoverage/S196_hexaploid_10Xcoverage.txt";

        /**
         * 提取LuLab JiaoLab 的10X数据，给康李鹏和治梁用
         */
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/fastcall/abd";
        String outfileDirS ="/data4/home/aoyue/vmap2/genotype/vcf_S306_10X_hexaploid_LuLabJiaoLab";
        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20201111";
        String taxaListS = "/data4/home/aoyue/vmap2/genotype/taxaList_toKang_withOnlyBamID_20201111.txt";



        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_hexaploid.vcf").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaListS + " > " + logfileS + " 2>&1 &");

//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.0.vcf.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.0_hexaploid_S196.vcf").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaListS + " > " + logfileS + " 2>&1 &");

            String infileS = new File(infileDirS,"chr" + chrArr[i] + ".ABDgenome.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_10X_fastcall_hexaploid.vcf").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            System.out.println("java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaListS + " > " + logfileS + " 2>&1 &");


        }

        //#/****************************** tetraploid ******************************/#
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1";
//        String outfileDirS ="/data4/home/aoyue/vmap2/genotype/mergedVCF/104_VCFbyPop/001_byPloidy/002_tetraploid";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20200526";
//        String taxaListS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/EmmerWheat_S187.txt";
//
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        for (int i = 0; i < chrArr.length; i++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_tetraploid.vcf").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaListS + " > " + logfileS + " 2>&1 &");
//        }

        //#/****************************** diploid ******************************/#
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1";
//        String outfileDirS ="/data4/home/aoyue/vmap2/genotype/mergedVCF/104_VCFbyPop/001_byPloidy/001_diploid/";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20200526";
//        String taxaListS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/Ae.tauschii_S36.txt";
//
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
//        for (int i = 0; i < chrArr.length; i++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_diploid.vcf").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaListS + " > " + logfileS + " 2>&1 &");
//        }

    }

    public void scriptFilterSNPtoBi(){
        //027_filterSNPtoBi_singleThread.jar
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_VMap2.0";
        String outfileDirS ="/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1";
        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log0605";
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};

//        for (int i = 0; i < chrArr.length; i++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.0.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1.vcf").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 027_filterSNPtoBi_singleThread.jar " + infileS + " " + outfileS + " > " + logfileS  + "" );
//        }

        SplitScript.splitScript2("/Users/Aoyue/Documents/sh.sh",7,12);

    }


    public void generateBiSNP(String infileS, String outfileS){

        File f = new File(infileS);

        try {
            String chr = f.getName().substring(3, 6); //提取染色体号 001
            BufferedReader br = AoFile.readFile(f.getAbsolutePath());
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            int biallelicNum = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                } else {
                    String kk = temp.substring(0,150);
                    String alt = PStringUtils.fastSplit(kk).get(4);
                    if (alt.contains(",")) continue;
                    if (alt.equals("D")) continue;
                    if (alt.equals("I")) continue;
                    biallelicNum++;
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(chr + "\t" + biallelicNum);
            System.out.println(f.getName() + " is completed at " + outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void calSite(){
        CountSites.countSites_fromVCF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subset/Dsubgenome.vcf.gz");

    }
    public void mergeTaxaDepth(){
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth_abd.summary.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth_ab.summary.txt";
        String infileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth_d.summary.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth.summary.txt";
        File f1 = new File(infileS1);
        File f2 = new File(infileS2);
        File f3 = new File(infileS3);
        File[] fs = {f1,f2,f3};
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            /**
             * 需要改动,header的名字可以自定义
             */
            //read header
            bw.write(br.readLine() + "\tGenomeType");
            bw.newLine();

            int cnttotal = 0;
            //read context
            int cnt = 0;
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                /**
                 * 需要改动
                 */
                String name = fs[i].getName();
                String group = name.split(".summary")[0].split("_")[1];
                if (group.equals("abd")){
                    group="AABBDD";
                }else if (group.equals("ab")){
                    group="AABB";
                }else if (group.equals("d")){
                    group="DD";
                }
                br = AoFile.readFile(infileS);
                br.readLine(); //read header
                String temp = null;

                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    cnttotal++;
                    String id = PStringUtils.getNDigitNumber(3,cnt);
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append("vmap2_").append(id).append("\t").append(l.get(2)).append("\t").append(group);
                    bw.write(sb.toString());
                    bw.newLine();
                }
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

    public void mkDepthSummary () {
        //思想：对每一个taxa做统计，每读进一个taxa，就统计所有深度的和，再除以表格的行数，也即是抽查的位点数，最终得出总得深度除以位点数，得出平均每个位点的深度。

//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/abd";
//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth_abd.summary.txt";

//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/ab";
//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth_ab.summary.txt";

        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/d";
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/taxaDepth_d.summary.txt";

//        String taxaDepthDirS = "/Users/Aoyue/Downloads/d";
//        String taxaSummaryFileS = "/Users/Aoyue/Downloads/d/taxaDepth_d.summary.txt";

        File[] fs = AoFile.getFileArrayInDir(taxaDepthDirS);
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(taxaSummaryFileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                RowTable t = new RowTable (fs[i].getAbsolutePath());
                String taxaName = String.valueOf(t.getHeader().get(0)).replaceFirst("_siteDepth", "");
                double value = 0;
                for (int j = 0; j < t.getRowNumber(); j++) {
                    value+=t.getCellAsDouble(j, 0);
                }
                double dd = (double)value/t.getRowNumber();
                bw.write(taxaName+"\t"+String.valueOf(i+1)+"\t"+String.format("%.4f", dd));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkDepthOfVMapII(){

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/ABDsubgenome_hexa.vcf.gz";
//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/abd";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/ABsubgenome_tetra.vcf.gz";
//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/ab";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/Dsubgenome_diploid.vcf.gz";
        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/004_taxaDepth/d";

//        String infileS = "/Users/Aoyue/Downloads/chr005_subset.vcf.gz";
//        String taxaDepthDirS = "/Users/Aoyue/Downloads/d";


        try {
            String temp = null;
            BufferedReader br = AoFile.readFile(infileS);
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l = PStringUtils.fastSplit(temp, "\t");
            String[] taxa = new String[l.size()-9]; //确定 taxa 的列表
            for (int i = 0; i < taxa.length; i++) {
                taxa[i] = l.get(i+9);
            }
            TIntArrayList[] depthList = new TIntArrayList[taxa.length]; //每个taxa建立一个深度集合
            for (int i = 0; i < taxa.length; i++) depthList[i] = new TIntArrayList(); //并初始化
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++; // 对snp开始计数
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt)+" lines");
                l = PStringUtils.fastSplit(temp, "\t");
                for (int i = 0; i < taxa.length; i++) {
                    String genoS = l.get(i+9);
                    if (genoS.startsWith(".")) {
                        depthList[i].add(0);
                        continue;
                    }
                    List<String> ll = PStringUtils.fastSplit(genoS, ":");
                    List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                    int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1));
                    depthList[i].add(depth);
                }
            }
            for (int i = 0; i < taxa.length; i++) {
                String outfileS = new File (taxaDepthDirS, PStringUtils.getNDigitNumber(3, i+1)+"depth.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(taxa[i]+"_siteDepth");
                bw.newLine();
                int[] depth = depthList[i].toArray();
                for (int j = 0; j < depth.length; j++) {
                    bw.write(String.valueOf(depth[j]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }



    public void mergeTxtandAddGroup() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/003_siteDepth";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/003_siteDepth/merge/site_depth.txt.gz";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/003_siteDepth/bySub";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/003_siteDepth/bySub/merge/site_depth_bySubgenome.txt.gz";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/004_siteDepth/001_byPloidy";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/004_siteDepth/002_mergeAndAddsub/site_depth_byPloidy.txt.gz";

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            /**
             * 需要改动,header的名字可以自定义
             */
            //read header
            bw.write(br.readLine() + "\tPloidy");
//            bw.write(br.readLine() + "\tSub");
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                /**
                 * 需要改动
                 */
                String name = fs[i].getName();
//                String group = name.split("subgenome")[0];
                String group = name.split("_depth")[0];


                br = AoFile.readFile(infileS);
                br.readLine(); //read header
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp).append("\t").append(group);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
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

    public void statVcfDepth_SD() {

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/003_siteDepth";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subset";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/003_siteDepth/bySub";

        // Indel的不同群体的深度和sd信息

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/000_VCFwithonlyIndel";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/004_siteDepth/001_byPloidy";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f  -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0] + "_depth.txt.gz").getAbsolutePath();

            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String[] taxa = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            try {
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##"))continue;
                    if(temp.startsWith("#C")){
                        List<String> linetaxa = PStringUtils.fastSplit(temp, "\t");
                        bw.write(linetaxa.get(0).replaceFirst("#", "") + "\t" + linetaxa.get(1) + "\t" + "AverageDepth\tSD");
                        bw.newLine();
                        taxa = new String[linetaxa.size() - 9];
                        continue;
                    }
                    cnt++; // 对snp开始计数
                    if (cnt % 1000000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    String chr = l.get(0);
                    String pos = l.get(1);
                    TDoubleArrayList depthList = new TDoubleArrayList();
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
                    sb.append(chr).append("\t").append(pos).append("\t").append(String.format("%.6f", relativeMean)).append("\t").append(String.format("%.6f", sd));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(infileS + " is calculated well done");

        });

    }

    public void getBinTable(){
        //*************** model *******************//
//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.05;
//        double max =0.5;
//        String outfileS = "";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/002_merge/001_site_QC.txt.gz";
//        AoFile.readheader(infileS);

//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.05;
//        double max = 0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/003_binTable/maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);

//        int indexGroup = 0;
//        int indexValue = 1;
//        double window = 0.05;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/003_binTable/site_heter.txt";
//        Bin.frequency2_byGroup(infileS,indexGroup,indexValue,window,window,outfileS);

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/002_merge/001_taxa_QC.txt.gz";
//        AoFile.readheader(infileS);
//        int indexGroup = 3;
//        int indexValue = 1;
//        double window = 0.01;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/003_binTable/taxa_heter.txt";
//        Bin.frequency2_byGroup(infileS,indexGroup,indexValue,window,window,outfileS);


        /**
         * Indel 的分Bin
         */
        //先检查header 0	GenomeType          //1	HeterozygousProportion          //2	MissingRate          //3	Maf
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/002_merge/001_site_QC.txt.gz";
//        AoFile.readheader(infileS);

        //maf的分bin
//        int indexGroup = 0;
//        int indexValue = 3;
//        double window = 0.009;
//        double max = 0.5;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/003_binTable/maf.txt";
//        Bin.frequency_byGroup(infileS,indexGroup,indexValue,max,window,window,outfileS);

        //heter的分bin
//        int indexGroup = 0;
//        int indexValue = 1;
//        double window = 0.05;
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/003_binTable/site_heter.txt";
//        Bin.frequency2_byGroup(infileS,indexGroup,indexValue,window,window,outfileS);


        //先检查 taxa质控的内容 0	Taxa        //1	HeterozygousProportion        //2	MissRate        //3	GenomeType
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/002_merge/001_taxa_QC.txt.gz";
        AoFile.readheader(infileS);
        int indexGroup = 3;
        int indexValue = 1;
        double window = 0.01;
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/003_binTable/taxa_heter.txt";
        Bin.frequency2_byGroup(infileS,indexGroup,indexValue,window,window,outfileS);


    }


    /**
     * VCF quality control
     * 先对合并的VCF进行群体的变异数目统计，再进行拆分群体 vcf
     * step 1: sample data all genome 0.001 即241K的数据 and then merge all
     * step 2:
     * step 2: abd ab d   Site: MAF miss heter depth  Taxa: Miss Heter
     */


    /**
     * 将生成的以site和taxa为单位的质控的结果进行合并，使六倍体，四倍体，二倍体在一个文件中
     */
    public void mergeCheckFile(){

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/002_merge";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/001_ori";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/007_IndelQC/002_merge";

        // change every time
//        String suffix = "_site_QC.txt.gz";
        String suffix = "_taxa_QC.txt.gz";

        String outfileS = new File(outfileDirS,"001" + suffix).getAbsolutePath();
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs,suffix);
        AoFile.mergeTxt_byFileArray(fs,outfileS);
    }

    public void QC(){
        String abdinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/ABDsubgenome_hexa.vcf.gz";
        String abinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/ABsubgenome_tetra.vcf.gz";
        String dinfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF/Dsubgenome_diploid.vcf.gz";

        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/001";

        this.checkQuality(abinfileS,outfileDirS,"AB");
        this.checkQuality(abdinfileS,outfileDirS,"ABD");
        this.checkQuality(dinfileS,outfileDirS,"D");


    }

    private void checkQuality (String infileS, String outDirS, String genomeType) {
        GenotypeGrid gt = new GenotypeGrid(infileS,GenoIOFormat.VCF_GZ);
        TDoubleArrayList missingSite = new TDoubleArrayList(); // calculation 1
        TDoubleArrayList hetSite = new TDoubleArrayList(); // calculation 2
        TDoubleArrayList maf = new TDoubleArrayList(); // calculation 3
        TDoubleArrayList missingTaxon = new TDoubleArrayList(); // calculation 4
        TDoubleArrayList hetTaxon = new TDoubleArrayList(); // calculation 5

        for (int i = 0; i < gt.getSiteNumber(); i++) {
            missingSite.add(((double) gt.getMissingNumberBySite(i)/gt.getTaxaNumber()));
            hetSite.add(gt.getHeterozygousProportionBySite(i));
            maf.add(gt.getMinorAlleleFrequency(i));

        }
        for (int i = 0; i < gt.getTaxaNumber(); i++) {
            missingTaxon.add((double)gt.getMissingNumberByTaxon(i)/gt.getSiteNumber());
            hetTaxon.add(gt.getHeterozygousProportionByTaxon(i));
        }

        String siteQCfileS = new File(outDirS,  genomeType + "_site_QC.txt.gz").getAbsolutePath();
        String taxaQCFileS = new File (outDirS, genomeType + "_taxa_QC.txt.gz").getAbsolutePath();

        try {
            BufferedWriter bw = AoFile.writeFile(siteQCfileS);
            bw.write("GenomeType\tHeterozygousProportion\tMissingRate\tMaf");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < missingSite.size(); i++) {
                sb.setLength(0);
                sb.append(genomeType).append("\t").append(hetSite.get(i)).append("\t").append(missingSite.get(i)).append("\t").append(maf.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }


        try {
            BufferedWriter bw = AoFile.writeFile(taxaQCFileS);
            bw.write("Taxa\tHeterozygousProportion\tMissRate\tGenomeType");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < gt.getTaxaNumber(); i++) {
                sb.setLength(0);
                sb.append(gt.getTaxonName(i)).append("\t").append(hetTaxon.get(i)).append("\t").append(missingTaxon.get(i)).append("\t").append(genomeType);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void mergeVCF(){

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001"; //这是做的测试！
//        String outfileDirS = "/Users/Aoyue/Documents/test";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/015_subesetVCF/singleChr";  //这是做的测试！
//        String infileDirS = "/Users/Aoyue/Documents/001_singleChr0.001"; //这是做的测试！
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF";

//        String infileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/001_fromVMap2.0_singleChr0.001";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/002_mergeVCFtoSub";

//        new CountSites().mergeVCFbysubgenome(infileDirS,outfileDirS);
//        new CountSites().mergeVCFtoAB_Dsubgenome(infileDirS,outfileDirS);
//        new CountSites().mergeVCFtoAandBsubgenome(infileDirS,outfileDirS);

//        CountSites.mergeVCFtoABsubgenome(infileDirS,outfileDirS);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/001_subsetVCF";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/temp/ABDsubgenome_hexa.vcf";
//        new CountSites().mergesubsetVCF(infileDirS,outfileS);
    }


    /**
     * Fix:VMapII
     * step 1: modify GL yo PL
     * step 2: rename VCF name like chr001_vmap2.0.vcf
     * step 3: bgzip vcf file and make index for vcf
     * step 4: make readme file for vmap2.0
     */

    public void sortTaxaName(){

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";


        String[] taxaArray = AoFile.getStringArraybyList_withoutHeader(infileS,0);
        Arrays.sort(taxaArray);
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            for (int i = 0; i < taxaArray.length; i++) {
                bw.write(taxaArray[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }


    public void bgzip() { // 重定向

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("bgzip chr" + chr + "_vmap2.0.vcf && tabix -p vcf chr" + chr + "_vmap2.0.vcf.gz &");
            System.out.println("bgzip -@ 4 chr" + chr + "_vmap2.0_hexaploid_S196.vcf && tabix -p vcf chr" + chr + "_vmap2.0_hexaploid_S196.vcf.gz &");
        }
    }

    public void modifyVMap2(){

        //sed -i '6s/GL/PL/' chr004_occu2_maf0.01_miss0.2.vcf
        // /pattern/s/pattern1/pattern2/： 上文的意思是只替换第6行的内容
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("sed -i '6s/GL/PL/' chr" + chr + "_occu2_maf0.01_miss0.2.vcf &");

            //chr042_occu2_maf0.01_miss0.2.vcf
            System.out.println("mv chr" +  chr + "_occu2_maf0.01_miss0.2.vcf chr" + chr + "_vmap2.0.vcf" );
        }
    }

    //chr041_vmap2.1.vcf.gz

    public void identifyFinalResult(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/003_occurenceCount/chr001_occu2_maf0.01_miss0.2.vcf";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
//            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                int cnt1 = Integer.parseInt(l.get(4));
                int cnt2 = Integer.parseInt(l.get(5));
                if(cnt1 <2 && cnt2 <2){
                    System.out.println(temp);
                }
            }
            br.close();

            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //验证都通过！！！！！！！
    }

    private String getSubgenomeInfoTestOccurrence(String[] PopGenoArray, String altList) {
//        int   dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
//        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
//            for (int j = 0; j < temList.size(); j++) {
//                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
//                dp += c; //dp是总深度
//                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
//            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
            int index1 = Integer.parseInt(temList.get(0)); //
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
            if (index1 != index2) {
                ht++;
            }
        }

        //计算 0/0 0/1 1/1 个数
        int cnt00 = gnCnt[0][0];
        int cnt01 = gnCnt[0][1] + gnCnt[1][0];
        int cnt11 = gnCnt[1][1];
        int cntmajor = cnt00 + cnt01;
        int cntminor = cnt01 + cnt11;
        int cntOccurrence = 0;
        if (cntmajor > cntminor){
            cntOccurrence = cntminor;
        }else{
            cntOccurrence = cntmajor;
        }

        //计算位点缺失率
//        nz = PopGenoArray.length - nz;
        float missRate = (float) ((double) nz/PopGenoArray.length);

        //计算 maf 和 aaf
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }
        float aaf = (float) ((double) acCnt[1] / sum);

        StringBuilder sb = new StringBuilder();
//        sb.append(String.format("%.4f", aaf)).append(",").append(String.format("%.4f", maf)).append(",").append(String.format("%.4f",missRate)); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D
        sb.append(String.format("%.4f", aaf)).append(",").append(String.format("%.4f", maf)).append(",").append(String.format("%.4f",missRate)).append(",").append(cntOccurrence); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D

        return sb.toString();
    }

    /**
     * 从上一步中获取的那些 不满足 maf >= 0.01 但是 minor 会出现在2个个体中的VCF，
     * 进行验证，打印出 每个群体 occur 出现次数
     */
    public void getOccurrenceByTaxa(){

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/002_occurrenceFile";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/003_occurenceCount";

        List<File> fList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fList);

        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";
        String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";
        String tetraFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";

        int occu = 2;
        float mafThresh = (float) 0.01;
        float missingThresh = (float) 0.2;



        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(hexaFileS,0);
        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(tetraFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(diFileS,0);
        Arrays.sort(abdTaxa);
        Arrays.sort(dTaxa);
        Arrays.sort(abTaxa);


//        File f = new File(infileS);

        fList.stream().forEach(f -> {
            StringBuilder s = new StringBuilder();
            s.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".vcf");
            String outfileS = new File(outfileDirS, s.toString()).getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);


            List<Integer> indexABD = new ArrayList<>();
            List<Integer> indexABorD = new ArrayList<>();
            String[] taxaABorDArray = null;
            String aaf = null;  //在注释文件中是写AAF_D 还是AAF_AB
            String annoHeader = null; //Header 中是 四倍体还是六倍体
            if (subgenome.equals("D")) {
                taxaABorDArray = dTaxa;
                aaf = "AAF_D";
                annoHeader = this.annotationHeader_Dsub();
            } else {
                taxaABorDArray = abTaxa;
                aaf = "AAF_AB";
                annoHeader = this.annotationHeader_ABsub();
            }

            System.out.println("Chr\tTotalSNP Num\tBiallelic Num\tTriallelic Num\tDeletion Num\tInsertion Num\tIndel Num");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                List<String> l = new ArrayList<>();
                int cntSNP = 0; //totalSNP
                int cntkept = 0;

                while ((temp = br.readLine()) != null) {
                    //***********************************************************//
                    if (temp.startsWith("##")) continue;

                    //***********************************************************//

                    if (temp.startsWith("#CHROM")) {
                        l = PStringUtils.fastSplit(temp);

                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index1 = Arrays.binarySearch(abdTaxa, taxon);
                            int index2 = Arrays.binarySearch(taxaABorDArray, taxon);

                            if (index1 > -1) {
                                indexABD.add(i);
                            }
                            if (index2 > -1) {
                                indexABorD.add(i);
                            }
                        }
                        Collections.sort(indexABD);
                        Collections.sort(indexABorD);
                    }
                    if (!temp.startsWith("#")) { //
                        cntSNP++;
                        l = PStringUtils.fastSplit(temp);
                        String altList = l.get(4);
                        int pos  = Integer.parseInt(l.get(1));
                        List<String> lgeno = new ArrayList<>();
                        List<String> lABDGeno = new ArrayList<>();
                        List<String> lABorDGeno = new ArrayList<>();


                        for (int i = 9; i < l.size(); i++) {
                            lgeno.add(l.get(i));
                        }
                        for (int i = 0; i < indexABD.size(); i++) {
                            lABDGeno.add(l.get(indexABD.get(i)));
                        }
                        for (int i = 0; i < indexABorD.size(); i++) {
                            lABorDGeno.add(l.get(indexABorD.get(i)));
                        }

                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                        String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                        String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                        String occur1 = this.getSubgenomeInfoTestOccurrence(hexaGenoArray, altList).split(",")[3];
                        String occur2 = this.getSubgenomeInfoTestOccurrence(ABorDGenoArray, altList).split(",")[3];

                        bw.write(l.get(0) + "\t" + l.get(1) + "\t" +  l.get(3)  + "\t" + l.get(4) + "\t" + occur1 + "\t" + occur2);
                        bw.newLine();

                    } //
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println( cntSNP + "\ttotal\t" + cntkept + " kept is completed at " + outfileS);

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });

    }


    /**
     * 检查哪些 MAF 小于 0.01 但是又属于 occurrence 大于 2 （即在2个个体中出现 minor allele 的位点）
      */
    public void getOccurrenceVCF(){

        String posMAFfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/chr001_occu2_maf0.01_miss0.2.txt";
        String posOccurrenceFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/occurrence/chr001_occu2_maf0.01_miss0.2.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/001_data/chr001.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/002_occurrenceFile/chr001.vcf";

        TIntArrayList posmaf = AoFile.getTIntList_withoutHeader(posMAFfileS,0);
        TIntArrayList posoccurr = AoFile.getTIntList_withoutHeader(posOccurrenceFileS,0);
        posmaf.sort();
        posoccurr.sort();

        try {

            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    int index = posmaf.binarySearch(pos);
                    int index2 = posoccurr.binarySearch(pos);
                    if (index > -1) continue;
                    if (index2 < 0) continue;
                    bw.write(temp);
                    bw.newLine();
                    cnt++;

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


    /**
     * 提取VCF文件的某几列信息
     */
    public void extractField(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_out2_byAomethod";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/004_extractColumn";
        String[] chrArr ={"001","005"};

        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_occu2_maf0.01_miss0.2.vcf").getAbsolutePath();
            String outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();
            int[] column = {1};
            AoFile.extractFileColumn(infileS,"#",column,outfileS);
        }
    }

    public void script(){

        //模板
//        String infileDirS = ""; //bcftools合并后的文件
//        String outfileDirS = ""; //进行maf miss occurrence 过滤的文件
//        String logDirS = ""; // 日志

//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF";
//        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_VMap2.0";
//        String logDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/log/004";

        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/106_MergedIndel_fromHapScanner";
        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/107_VMap2.0_Indel";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/029_hapScanner_onlyIndel/log_filterMafMissOccurrence/";

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};

//        String[] chrArr = {"010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String[] chrArr = {"003","004","005","006","007","008","009"};

        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + ".vcf").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(infileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -Xms200g -Xmx200g -jar 044_filterMafOccurrenceMiss.jar " + infileS + " " + outfileDirS + " > " + logfileS + " 2 >&1 &" );
            System.out.println("java -jar 044_filterMafOccurrenceMiss.jar " + infileS + " " + outfileDirS + " > " + logfileS + " 2 >&1 &" );

        }
    }

    public void filter_singleThread (String infileS, String outfileDirS) {

//        String outputVCFDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_MAF0.01";
//        String outputVCFDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/004_filterMAFmissOccurrence";


        String ABTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/EmmerWheat_S187.txt";
        String ABDTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/BreadWheat_S420.txt";
        String DTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/Ae.tauschii_S36.txt";

        int occu = 2;
        float mafThresh = (float)0.01;
        float missingThresh = (float)0.2;
        File f = new File(infileS);

        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(ABTaxaFileS,0);
        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(ABDTaxaFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(DTaxaFileS,0);

        StringBuilder sb = new StringBuilder();
        sb.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".vcf");
        String outfileS = new File (outfileDirS, sb.toString()).getAbsolutePath();
        int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
        String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);

        //总的VCF文件，获取 gt 和 allTaxa
        GenotypeGrid gt = new GenotypeGrid(f.getAbsolutePath(), GenoIOFormat.VCF);
        gt.sortByTaxa();
        String[] allTaxa = gt.getTaxaNames();
        /**
         * 建立 2 个群体的index
         */
        // 将VCF文件中的两个群体拆开
        int[][] taxaIndices = new int[2][];
        String[][] subTaxa = new String[2][];
        GenotypeGrid[] gts = new GenotypeGrid[2]; //群体的 genotypeTable
        if (subgenome.equals("D")) {
            taxaIndices[0] = new int[abdTaxa.length];
            taxaIndices[1] = new int[dTaxa.length];
            subTaxa[0] = abdTaxa;
            subTaxa[1] = dTaxa;
        }
        else {
            taxaIndices[0] = new int[abTaxa.length];
            taxaIndices[1] = new int[abdTaxa.length];
            subTaxa[0] = abTaxa;
            subTaxa[1] = abdTaxa;
        }

        for (int i = 0; i < taxaIndices.length; i++) { // 第一维，确定 TaxaIndices是哪个群体的
            for (int j = 0; j < taxaIndices[i].length; j++) { //第二维，确定群体内每个 taxa 在VCF文件中的 allTaxa中的索引
                taxaIndices[i][j] = Arrays.binarySearch(allTaxa, subTaxa[i][j]);
                if (taxaIndices[i][j] < 0) System.out.println(subTaxa[i][j]);
            }
        }

        // 从上文建立的 taxaIndex里，选出 单个群体，new gts[i]
        for (int i = 0; i < gts.length; i++) {
            gts[i] = GenotypeOperation.getSubsetGenotypeByTaxon(gt, taxaIndices[i]);
        }
        TIntArrayList posList = new TIntArrayList(); //这是我们将要保留的sites
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            double[] missing = new double[2];
            for (int j = 0; j < 2; j++) { //第一步，获取2个群体的 missing rate
                missing[j] = (double)gts[j].getMissingNumberBySite(i)/gts[j].getTaxaNumber();
            }
            if (missing[0] > missingThresh && missing[1] > missingThresh) continue; //如果缺失率都大于0.2的话，该位点不保留。

            for (int j = 0; j < gts.length; j++) { //判断2个群体的maf值和occurrence 值

                double maf = gts[j].getMinorAlleleFrequency(i);
                if (!(maf < mafThresh)) { //
                    posList.add(gt.getPosition(i));
                    break;
                }
                else if (!(gts[j].getAlternativeAlleleOccurrenceBySite(i) < occu)) {
                    posList.add(gt.getPosition(i));
                    break;
                }
            }
        }
        int[] positions = posList.toArray(); //*********** 该pos库用作下文搜索用   *****************//
        Arrays.sort(positions);
        System.out.println(positions.length + " sites were kept in " + new File(outfileS).getName());
        System.out.println("***********************************************************************");
        System.out.println("***********************************************************************");


        /**
         * begin to write the kept site in vcf file
         */
        try{
            List<Integer> indexABD = new ArrayList<>();
            List<Integer> indexABorD = new ArrayList<>();
            String[] taxaABorDArray = null;
            String aaf = null;  //在注释文件中是写AAF_D 还是AAF_AB
            String annoHeader = null; //Header 中是 四倍体还是六倍体
            if (subgenome.equals("D")){
                taxaABorDArray = dTaxa;
                aaf = "AAF_D";
                annoHeader = this.annotationHeader_Dsub();
            }else{
                taxaABorDArray = abTaxa;
                aaf = "AAF_AB";
                annoHeader = this.annotationHeader_ABsub();
            }

            BufferedReader br = AoFile.readFile(f.getAbsolutePath());
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(annoHeader); bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            while((temp=br.readLine()) != null){
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#CHROM")){
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(abdTaxa, taxon);
                        int index2 = Arrays.binarySearch(taxaABorDArray, taxon);
                        if (index1 > -1) {
                            indexABD.add(i);
                        }
                        if (index2 > -1) {
                            indexABorD.add(i);
                        }
                    }
                    Collections.sort(indexABD);
                    Collections.sort(indexABorD);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    int kk = Arrays.binarySearch(positions,pos);
                    if (kk < 0) continue;
                    String altList = l.get(4);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lABDGeno = new ArrayList<>();
                    List<String> lABorDGeno = new ArrayList<>();

                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexABD.size(); i++) {
                        lABDGeno.add(l.get(indexABD.get(i)));
                    }
                    for (int i = 0; i < indexABorD.size(); i++) {
                        lABorDGeno.add(l.get(indexABorD.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                    String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                    String INFO = this.getInfo(genoArray, altList);
                    String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
                    String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];


                    StringBuilder sbb = new StringBuilder();
                    for (int i = 0; i < 7; i++) {
                        sbb.append(l.get(i)).append("\t");
                    }

                    sbb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";").append(aaf).append("=").append(ABorDAAF).append("\tGT:AD:GL");
                    for (int i = 9; i < l.size(); i++) {
                        sbb.append("\t").append(l.get(i));
                    }
                    bw.write(sbb.toString());
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(outfileS + " is completed.");

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public String annotationHeader_Dsub(){
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);

        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n" +
                "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
                "##fileDate=" + S.split(" ")[0] + "\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n" +
                "##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n" +
                "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n" +
                "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n" +
                "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n" +
                "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n" +
                "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n" +
                "##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">\n" +
                "##INFO=<ID=AAF_D,Number=1,Type=Float,Description=\"Alternative allele frequency on diploid Aegilops tauschii\">\n" +  //这一行根据 D 还是 AB 有变化
                "##ALT=<ID=D,Description=\"Deletion\">\n" +
                "##ALT=<ID=I,Description=\"Insertion\">\n" +
                "##Species=Wheat\n" +
                "##ReferenceGenome=iwgsc_refseqv1.0\n" +
                "##VariantsMapVersion=\"vmap2\"");

        return sb.toString();

    }

    public String annotationHeader_ABsub(){
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);

        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n" +
                "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
                "##fileDate=" + S.split(" ")[0] + "\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n" +
                "##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n" +
                "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n" +
                "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n" +
                "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n" +
                "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n" +
                "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n" +
                "##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">\n" +
                "##INFO=<ID=AAF_AB,Number=1,Type=Float,Description=\"Alternative allele frequency on tetraploid emmer wheat\">\n" +   //这一行根据 D 还是 AB 有变化
                "##ALT=<ID=D,Description=\"Deletion\">\n" +
                "##ALT=<ID=I,Description=\"Insertion\">\n" +
                "##Species=Wheat\n" +
                "##ReferenceGenome=iwgsc_refseqv1.0\n" +
                "##VariantsMapVersion=\"vmap2\"");

        return sb.toString();
    }


    public void filter_parallel () {

//        String inputVCFDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF";
//        String outputVCFDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_MAF0.01";
//        String ABTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/EmmerWheat_S187.txt";
//        String ABDTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/BreadWheat_S420.txt";
//        String DTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/Ae.tauschii_S36.txt";


        String inputVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/001_data";
        String outputVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out";
//        String outputVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/occurrence";

        String ABTaxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";
        String ABDTaxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";
        String DTaxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";

        int occu = 2;
        float mafThresh = (float)0.01;
        float missingThresh = (float)0.2;
        List<File> fList = AoFile.getFileListInDir(inputVCFDirS);
        Collections.sort(fList);

        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(ABTaxaFileS,0);
        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(ABDTaxaFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(DTaxaFileS,0);

        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".txt");
            String outfileS = new File (outputVCFDirS, sb.toString()).getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);

            //总的VCF文件，获取 gt 和 allTaxa
            GenotypeGrid gt = new GenotypeGrid(f.getAbsolutePath(), GenoIOFormat.VCF_GZ);
            gt.sortByTaxa();
            String[] allTaxa = gt.getTaxaNames();
            /**
             * 建立 2 个群体的index
             */
            // 将VCF文件中的两个群体拆开
            int[][] taxaIndices = new int[2][];
            String[][] subTaxa = new String[2][];
            GenotypeGrid[] gts = new GenotypeGrid[2]; //群体的 genotypeTable
            if (subgenome.equals("D")) {
                taxaIndices[0] = new int[abdTaxa.length];
                taxaIndices[1] = new int[dTaxa.length];
                subTaxa[0] = abdTaxa;
                subTaxa[1] = dTaxa;
            }
            else {
                taxaIndices[0] = new int[abTaxa.length];
                taxaIndices[1] = new int[abdTaxa.length];
                subTaxa[0] = abTaxa;
                subTaxa[1] = abdTaxa;
            }

            for (int i = 0; i < taxaIndices.length; i++) { // 第一维，确定 TaxaIndices是哪个群体的
                for (int j = 0; j < taxaIndices[i].length; j++) { //第二维，确定群体内每个 taxa 在VCF文件中的 allTaxa中的索引
                    taxaIndices[i][j] = Arrays.binarySearch(allTaxa, subTaxa[i][j]);
                    if (taxaIndices[i][j] < 0) System.out.println(subTaxa[i][j]);
                }
            }

            // 从上文建立的 taxaIndex里，选出 单个群体，new gts[i]
            for (int i = 0; i < gts.length; i++) {
                gts[i] = GenotypeOperation.getSubsetGenotypeByTaxon(gt, taxaIndices[i]);
            }
            TIntArrayList posList = new TIntArrayList(); //这是我们将要保留的sites
            for (int i = 0; i < gt.getSiteNumber(); i++) {
                double[] missing = new double[2];
                for (int j = 0; j < 2; j++) { //第一步，获取2个群体的 missing rate
                    missing[j] = (double)gts[j].getMissingNumberBySite(i)/gts[j].getTaxaNumber();
                }
                if (missing[0] > missingThresh && missing[1] > missingThresh) continue; //如果缺失率都大于0.2的话，该位点不保留。

                for (int j = 0; j < gts.length; j++) { //判断2个群体的maf值和occurrence 值

                    double maf = gts[j].getMinorAlleleFrequency(i);
                    if (!(maf < mafThresh)) { //
                        posList.add(gt.getPosition(i));
                        break;
                    }
                    else if (!(gts[j].getAlternativeAlleleOccurrenceBySite(i) < occu)) {
                        posList.add(gt.getPosition(i));
                        break;
                    }

//                    if (!(gts[j].getAlternativeAlleleOccurrenceBySite(i) < occu)) {
//                        posList.add(gt.getPosition(i));
//                        break;
//                    }

                }

            }
            int[] positions = posList.toArray();
            Arrays.sort(positions);
            System.out.println(positions.length + " sites were kept in " + new File(outfileS).getName());
            System.out.println("***********************************************************************");
            System.out.println("***********************************************************************");

            /**
             * 测试老师的MAF和我的MAF之间的差别
             *
             */
//            try{
//                BufferedWriter bw = AoFile.writeFile(outfileS);
//                for (int i = 0; i < positions.length; i++) {
//                    bw.write(String.valueOf(positions[i]));
//                    bw.newLine();
//                }
//                bw.flush();
//                bw.close();
//
//            }catch (Exception e) {
//                e.printStackTrace();
//                System.exit(1);
//            }



            /**
             * begin to write the kept site in vcf file
             */
            try{
                List<Integer> indexABD = new ArrayList<>();
                List<Integer> indexABorD = new ArrayList<>();
                String[] taxaABorDArray = null;
                String aaf = null;
                String annoHeader = null;
                if (subgenome.equals("D")){
                    taxaABorDArray = dTaxa;
                    aaf = "AAF_D";
                    annoHeader = this.annotationHeader_Dsub();
                }else{
                    taxaABorDArray = abTaxa;
                    aaf = "AAF_AB";
                    annoHeader = this.annotationHeader_ABsub();
                }

                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(annoHeader); bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while((temp=br.readLine()) != null){
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#CHROM")){
                        l = PStringUtils.fastSplit(temp);
                        bw.write(temp);
                        bw.newLine();
                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index1 = Arrays.binarySearch(abdTaxa, taxon);
                            int index2 = Arrays.binarySearch(taxaABorDArray, taxon);
                            if (index1 > -1) {
                                indexABD.add(i);
                            }
                            if (index2 > -1) {
                                indexABorD.add(i);
                            }
                        }
                        Collections.sort(indexABD);
                        Collections.sort(indexABorD);
                    }
                    if (!temp.startsWith("#")) {
                        l = PStringUtils.fastSplit(temp);
                        int pos = Integer.parseInt(l.get(1));
                        int kk = Arrays.binarySearch(positions,pos);
                        if (kk < 0) continue;
                        String altList = l.get(4);
                        List<String> lgeno = new ArrayList<>();
                        List<String> lABDGeno = new ArrayList<>();
                        List<String> lABorDGeno = new ArrayList<>();

                        for (int i = 9; i < l.size(); i++) {
                            lgeno.add(l.get(i));
                        }
                        for (int i = 0; i < indexABD.size(); i++) {
                            lABDGeno.add(l.get(indexABD.get(i)));
                        }
                        for (int i = 0; i < indexABorD.size(); i++) {
                            lABorDGeno.add(l.get(indexABorD.get(i)));
                        }

                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                        String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                        String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                        String INFO = this.getInfo(genoArray, altList);
                        String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
                        String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];


                        StringBuilder sbb = new StringBuilder();
                        for (int i = 0; i < 7; i++) {
                            sbb.append(l.get(i)).append("\t");
                        }

                        sbb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";").append(aaf).append("=").append(ABorDAAF).append("\tGT:AD:GL");
                        for (int i = 9; i < l.size(); i++) {
                            sbb.append("\t").append(l.get(i));
                        }
                        bw.write(sbb.toString());
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(outfileS + " is completed.");

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });

    }



//String infileS, String outfileDirS

    public void filterMafbyPopHTD() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/001_data";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_out2_byAomethod";

        List<File> fList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fList);

        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";
        String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";
        String tetraFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";

        int occu = 2;
        float mafThresh = (float) 0.01;
        float missingThresh = (float) 0.2;


        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(hexaFileS,0);
        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(tetraFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(diFileS,0);
        Arrays.sort(abdTaxa);
        Arrays.sort(dTaxa);
        Arrays.sort(abTaxa);


//        File f = new File(infileS);

        fList.stream().forEach(f -> {
            StringBuilder s = new StringBuilder();
            s.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".vcf");
            String outfileS = new File(outfileDirS, s.toString()).getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);


            List<Integer> indexABD = new ArrayList<>();
            List<Integer> indexABorD = new ArrayList<>();
            String[] taxaABorDArray = null;
            String aaf = null;  //在注释文件中是写AAF_D 还是AAF_AB
            String annoHeader = null; //Header 中是 四倍体还是六倍体
            if (subgenome.equals("D")) {
                taxaABorDArray = dTaxa;
                aaf = "AAF_D";
                annoHeader = this.annotationHeader_Dsub();
            } else {
                taxaABorDArray = abTaxa;
                aaf = "AAF_AB";
                annoHeader = this.annotationHeader_ABsub();
            }

            System.out.println("Chr\tTotalSNP Num\tBiallelic Num\tTriallelic Num\tDeletion Num\tInsertion Num\tIndel Num");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(annoHeader);
                bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cntSNP = 0; //totalSNP
                int cntkept = 0;

                while ((temp = br.readLine()) != null) {
                    //*********************** # section ************************************//
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#CHROM")) {
                        l = PStringUtils.fastSplit(temp);
                        bw.write(temp);
                        bw.newLine();

                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index1 = Arrays.binarySearch(abdTaxa, taxon);
                            int index2 = Arrays.binarySearch(taxaABorDArray, taxon);

                            if (index1 > -1) {
                                indexABD.add(i);
                            }
                            if (index2 > -1) {
                                indexABorD.add(i);
                            }
                        }
                        Collections.sort(indexABD);
                        Collections.sort(indexABorD);
                    }
                    //*********************** pos section ************************************//
                    if (!temp.startsWith("#")) { //
                        cntSNP++;
                        l = PStringUtils.fastSplit(temp);
                        String altList = l.get(4);
                        int pos  = Integer.parseInt(l.get(1));
                        List<String> lgeno = new ArrayList<>();
                        List<String> lABDGeno = new ArrayList<>();
                        List<String> lABorDGeno = new ArrayList<>();

                        for (int i = 9; i < l.size(); i++) {
                            lgeno.add(l.get(i));
                        }
                        for (int i = 0; i < indexABD.size(); i++) {
                            lABDGeno.add(l.get(indexABD.get(i)));
                        }
                        for (int i = 0; i < indexABorD.size(); i++) {
                            lABorDGeno.add(l.get(indexABorD.get(i)));
                        }

                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                        String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                        String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                        String hexaMAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[1];
                        String diMAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[1];

                        //开始进行,maf判断
                        double hexamaf = Double.parseDouble(hexaMAF);
                        double ABorDmaf = Double.parseDouble(diMAF);

                        if (hexamaf >= 0.01 || (ABorDmaf >= 0.01)) {
                            cntkept++;

                            String INFO = this.getInfo(genoArray, altList);
                            String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
                            String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];

                            StringBuilder sb = new StringBuilder();
                            for (int i = 0; i < 7; i++) {
                                sb.append(l.get(i)).append("\t");
                            }
                            sb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";").append(aaf).append("=").append(ABorDAAF).append("\tGT:AD:GL");
                            for (int i = 9; i < l.size(); i++) {
                                sb.append("\t").append(l.get(i));
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println( cntSNP + "\ttotal\t" + cntkept + " kept is completed at " + outfileS);

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });

    }

    private String getInfo(String[] genoArray, String altList) {
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
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
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

    private String getSubgenomeInfo(String[] PopGenoArray, String altList) {
//        int   dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
//        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
//        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
//            for (int j = 0; j < temList.size(); j++) {
//                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
//                dp += c; //dp是总深度
//                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
//            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
//            int index1 = Integer.parseInt(temList.get(0)); //
//            int index2 = Integer.parseInt(temList.get(1));
//            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
//            if (index1 != index2) {
//                ht++;
//            }
        }
//        nz = PopGenoArray.length - nz;
        float missRate = (float) ((double) nz/PopGenoArray.length);

        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }
        float aaf = (float) ((double) acCnt[1] / sum);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%.4f", aaf)).append(",").append(String.format("%.4f", maf)).append(",").append(String.format("%.4f",missRate)); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D
        return sb.toString();
    }
}
