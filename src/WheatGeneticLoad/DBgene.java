package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.AoString;
import AoUtils.CalVCF;
import AoUtils.Triads.Standardization;
import AoUtils.Triads.Triadsgenes;
import daxing.common.RowTableTool;


import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;
import pgl.infra.anno.gene.GeneFeature;

import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.function.Predicate;

import pgl.infra.range.Range;

public class DBgene {

    public DBgene(){
//        this.getHexaploidAnnotation();
//        this.getTranscriptSum(); //多线程运行造成内存不足，故使用单个染色体进行运行
//        this.getTranscriptSum_bychr();
//        this.script_getTranscriptSum();
//        this.mergeTxt();

        /**
         *  开始处理cultivar 和 landrace 的
         */
//        this.getsubspeciesSNPAnnotation();
//        this.script_getTranscriptSum();
//        this.mergeTxt();
//        this.addGroupforgeneSummary();

        /**
         * TriadsID
         */
//        this.mkSpreadFormat_hexaploid();
//        this.getsubspeciesSNPAnnotation_Tetraploid_diploid();
//        this.script_getTranscriptSum();
//        this.mergeTxt();
//        this.mkSpreadFormat_tetraploid_diploid();
//        this.mergeSpreadTable_HexaploidTetraploidDiploid();
//        this.mergeSpreadTable();

        /**
         * 2020-11-19 update my data
         */
//        this.getHexaploidAnnotation_();
        this.script_getTranscriptSum();
//        this.mergeTxt();
//        this.mkSpreadFormat_hexaploid();
//        this.addTriadIDforEpigenomicMap();

        /**
         * 2021-04-15 update tetraploid and diploid
         */

//        this.getPseudoHexaploidAnnotation_(); //
//        this.script_getTranscriptSum();
//        this.mergeTxt();
//        this.mkSpreadFormat_hexaploid();








    }


    /**
     * 向表观图谱加入 triadsID,统一编号
     */
    public void addTriadIDforEpigenomicMap(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/source/Epigenomic_map_inWheat.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/004_epigenomic/001_EpigenomicMap_addTriadsID.txt.gz";

        Triadsgenes tg = new Triadsgenes();
        System.out.println("cnt" + tg.getTriadNum());
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = null;
            List<String> headList = PStringUtils.fastSplit(br.readLine());
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < headList.size(); i++) {
                sb.append(headList.get(i)).append("\t");
                i++;
            }
            header = sb.deleteCharAt(sb.length()-1).toString();
            bw.write("TriadID" + "\t");
            bw.write(header);bw.newLine();
            br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            String geneA;
            String geneA_02;

            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(temp);

                //判断8组geneid是否相等
//                for (int i = 0; i < 14; i++) {
//                    if (!l.get(i).equals(l.get(i+2))){
//                        System.out.println(l.get(i) + "\t" + l.get(i+2));
//                    }
//                    i++;
//                }
                geneA = l.get(0).split(";")[0]; //该基因为01版本的基因
                geneA_02 = AoString.getv11geneName(geneA);
                String triadID = tg.getTriadID(geneA_02);
                if(triadID == null || triadID == ""){  //!!!!! if there is no value, we should set the value as "NA".
                    triadID = "NA";
                }
                else{
                    sb.append(triadID);
                    for (int i = 0; i < l.size(); i++) {
                        i++;
                        sb.append("\t").append(l.get(i));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println( "epigenomic " + " count " + cnt);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 2021-04-15 周四 update
     */
    public void getPseudoHexaploidAnnotation_(){
//        String inputDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/";
//        String outDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/001_exonSNPAnnotation_hexaploid";

        String inputDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";
        String outDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/101_exonSNPAnnotation_pseudoHexaploid";

        List<File> files= AoFile.getFileListInDir(inputDir);
        String[] outNames=files.stream().map(File::getName).map(s->s.replaceAll("_anno.txt.gz", "_pseudohexaploid_anno.txt")).toArray(String[]::new);
        RowTableTool<String> rowTable;
        Predicate<List<String>> removed= l->Double.parseDouble(l.get(9))==0; //即AAF_ABD=0，过滤没有分离的位点
        for (int i = 0; i < files.size(); i++) {
            rowTable=new RowTableTool<>(files.get(i).getAbsolutePath());
            rowTable.removeIf(removed);
            rowTable.write(new File(outDir, outNames[i]).getAbsolutePath());
        }
    }

    /**
     * 2020-11-19 周四 update
     */
    public void getHexaploidAnnotation_(){
//        String inputDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/";
//        String outDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/001_exonSNPAnnotation_hexaploid";

        String inputDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";
        String outDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/001_exonSNPAnnotation_hexaploid";

        List<File> files= AoFile.getFileListInDir(inputDir);
        String[] outNames=files.stream().map(File::getName).map(s->s.replaceAll("_anno.txt.gz", "_hexaploid_anno.txt")).toArray(String[]::new);
        RowTableTool<String> rowTable;
        Predicate<List<String>> removed= l->Double.parseDouble(l.get(8))==0; //即AAF_ABD=0，过滤没有分离的位点
        for (int i = 0; i < files.size(); i++) {
            rowTable=new RowTableTool<>(files.get(i).getAbsolutePath());
            rowTable.removeIf(removed);
            rowTable.write(new File(outDir, outNames[i]).getAbsolutePath());
        }
    }


    /**
     *  求出3个表格的交集,即 landrace cultivar and Tetraploid Diploid
     */
    public void mergeSpreadTable(){
        //landrace
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/001_hexaploid/001_spreadTable_landraceCultivar/001_triadsLoad_landrace.txt";
        //cultivar
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/001_hexaploid/001_spreadTable_landraceCultivar/001_triadsLoad_cultivar.txt";
        // merged tetraploid diploid
        String infileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/003_spreadTable/001_triadsLoad_tetraploid_diploid.txt";

        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/004_alluvia/002_spreadTable_merged_alluvia_LR_CL.txt";
        List<String> l1 = new ArrayList<>();
        List<String> l2 = new ArrayList<>();
        List<String> l3 = new ArrayList<>();
        HashMap<String,String> hm1 = new HashMap<>();
        HashMap<String,String> hm2 = new HashMap<>();
        HashMap<String,String> hm3 = new HashMap<>();
        AoFile.readheader(infileS1);
        RowTable<String> t = new RowTable<>(infileS1);
        //1.kept the triad which are not the M000 model
        for (int i = 0; i < t.getRowNumber(); i++) {
            String triad = t.getCell(i,0);
            String region = t.getCell(i,4);
            if (region.equals("M000"))continue;
            hm1.put(triad,region);
            l1.add(triad);
        }
        System.out.println(l1.size() + " l1 size");

        t = new RowTable<>(infileS2);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String triad = t.getCell(i,0);
            String region = t.getCell(i,4);
            if (region.equals("M000"))continue;
            hm2.put(triad,region);
            l2.add(triad);
        }
        System.out.println(l2.size() + " l2 size");

        t = new RowTable<>(infileS3);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String triad = t.getCell(i,0);
            String region = t.getCell(i,4);
            if (region.equals("M000"))continue;
            hm3.put(triad,region);
            l3.add(triad);
        }
        System.out.println(l3.size() + " l3 size");
        l1.retainAll(l2);
        l1.retainAll(l3);
        System.out.println(l1.size() + " l1 and l2 and l3 intersection.");
        int a=3;

        try{
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("TriadID\tnonsynVSsynRatio_RegionOnAABB_DD\tnonsynVSsynRatio_RegionOnLandrace\tnonsynVSsynRatio_RegionOnCultivar\tFreq");
            bw.newLine();
            for (int i = 0; i < l1.size(); i++) {
                String triad = l1.get(i);
                String region1 = hm1.get(triad);
                String region2 = hm2.get(triad);
                String region3 = hm3.get(triad);
                bw.write(triad + "\t" + region3 + "\t" + region1+ "\t" + region2 + "\t1");
                bw.newLine();
            }
            bw.flush();
            bw.close();
            //3364 l1 size
            //4055 l2 size
            //9852 l3 size
            //2868 l1 and l2 and l3 intersection.
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     *  求出2个表格的交集
     */
    public void mergeSpreadTable_HexaploidTetraploidDiploid(){
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/001_hexaploid/test.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/003_spreadTable/001_triadsLoad_tetraploid_diploid.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/004_alluvia/001_spreadTable_merged_alluvia.txt";
        List<String> l1 = new ArrayList<>();
        List<String> l2 = new ArrayList<>();
        HashMap<String,String> hm1 = new HashMap<>();
        HashMap<String,String> hm2 = new HashMap<>();
        AoFile.readheader(infileS1);
        RowTable<String> t = new RowTable<>(infileS1);
        //1.kept the triad which are not the M000 model
        for (int i = 0; i < t.getRowNumber(); i++) {
            String triad = t.getCell(i,0);
            String region = t.getCell(i,4);
            if (region.equals("M000"))continue;
            hm1.put(triad,region);
            l1.add(triad);
        }
        System.out.println(l1.size() + " l1 size");
        t = new RowTable<>(infileS2);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String triad = t.getCell(i,0);
            String region = t.getCell(i,4);
            if (region.equals("M000"))continue;
            hm2.put(triad,region);
            l2.add(triad);
        }
        System.out.println(l2.size() + " l2 size");
        l1.retainAll(l2);
        System.out.println(l1.size() + " l1 and l2 intersection.");
        int a=3;

        try{
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("TriadID\tnonsynVSsynRatio_RegionOnAABBDD\tnonsynVSsynRatio_RegionOnAABB_DD\tFreq");
            bw.newLine();
            for (int i = 0; i < l1.size(); i++) {
                String triad = l1.get(i);
                String region1 = hm1.get(triad);
                String region2 = hm2.get(triad);
                bw.write(triad + "\t" + region1 + "\t" + region2 + "\t1");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 从最原始的geneSummary开始进行转换
     *
     */
    public void mkSpreadFormat_tetraploid_diploid(){
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/003_spreadTable/001_triadsLoad_tetraploid_diploid.txt";
        GeneDB genedb = new GeneDB();
        Triadsgenes tg = new Triadsgenes();
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("TriadID\tNonVsSynRatioA\tNonVsSynRatioB\tNonVsSynRatioD\tNonVsSynRatioRegion");
            bw.newLine();
            int cnt=0;
            int cntremaining =0;
            for (int i = 0; i < tg.getTriadNum(); i++) {
                cnt++;
                String triadID = tg.triadsList.get(i);
                String genea = tg.getGeneinAsub(triadID);
                String geneb = tg.getGeneinBsub(triadID);
                String gened = tg.getGeneinDsub(triadID);
                String ratioA = genedb.getNonVsSynRatio(genea);
                String ratioB = genedb.getNonVsSynRatio(geneb);
                String ratioD = genedb.getNonVsSynRatio(gened);
                //filter NA
                if (ratioA.startsWith("N") || ratioB.startsWith("N") || ratioD.startsWith("N")) continue;
                double[] ratiodABD = {Double.parseDouble(ratioA),Double.parseDouble(ratioB),Double.parseDouble(ratioD)};
                String region = Standardization.getNearestPointIndex(ratiodABD).getRegion();
                bw.write(triadID+"\t"+ratioA+"\t"+ratioB+"\t"+ratioD+"\t"+region);
                bw.newLine();
                cntremaining++;
            }
            System.out.println(cnt + " triads totally, " + cntremaining + " triads kept");
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据亚群的个体，结合merged VCF 文件，提取亚群的 snp annotation info
     *
     */
    public void getsubspeciesSNPAnnotation_Tetraploid_diploid(){
// need modify need modify need modify need modify
        //        String taxaList = "";
//        String outfileDirS = "";

        String vcfFileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony"; //chr001_SNP_anno.txt.gz

        //tetraploid
//        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_byPloidy/Tetraploid.txt";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/001_annotation/001";

        //diploid
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_byPloidy/Ae.tauschii.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/001_annotation/001";


        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        AoFile.readheader(fsList.get(0).getAbsolutePath());
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = "chr" + chrArr[i]; // = chr001
            String exonVCF = new File(vcfFileDirS,chr + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            //*********************************** step1: 获取群体的 pos 信息 ***********************************//
            TIntArrayList posl = CalVCF.extractVCFchrPos(exonVCF, taxaList);

            //******* step2: 根据annotation 库， pos列， pols目标集合， 输出文件获取 注释文件的子集 ******//
            String infileS = new File(infileDirS, chr + "_SNP_anno.txt.gz").getAbsolutePath(); //chr001_SNP_anno.txt.gz
            //需要修改，输出文件的名字 //需要修改，输出文件的名字 //需要修改，输出文件的名字
//            String outfileS = new File(outfileDirS, chr + "_SNP_anno_tetraploid.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS, chr + "_SNP_anno_diploid.txt.gz").getAbsolutePath();

            File out = AoFile.filterTxtLines(infileS,2,posl,outfileS);
        }
    }


    /**
     * 从最原始的geneSummary开始进行转换
     *
     */
    public void mkSpreadFormat_hexaploid(){
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/001_hexaploid/test.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/001_hexaploid/001_spreadTable_landraceCultivar/001_triadsLoad_landrace.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/001_hexaploid/001_spreadTable_landraceCultivar/001_triadsLoad_cultivar.txt";

        /**
         * 2020-11-19 update result
         */
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/003_traids_load/001_triadsLoad_nonsynVSsyn_hexaploid.txt"; //nonsyn vs syn
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/003_traids_load/002_triadsLoad_delVSsyn_hexaploid.txt"; // del VS syn
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/003_traids_load/003_triadsLoad_delFre_hexaploid.txt"; // del frequency

        /**
         * 2021-04-15 update result
         */

//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/103_traids_load/001_triadsLoad_nonsynVSsyn_pseudohexaploid.txt"; //nonsyn vs syn
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/103_traids_load/002_triadsLoad_delVSsyn_pseudohexaploid.txt"; // del VS syn
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/103_traids_load/003_triadsLoad_delFre_pseudohexaploid.txt"; // del frequency

        GeneDB genedb = new GeneDB(); //需要修改, 每批数据都有一个基因的summary，这里要更新
        Triadsgenes tg = new Triadsgenes();  //需要修改； 有3个版本
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            bw.write("TriadID\tNonVsSynRatioA\tNonVsSynRatioB\tNonVsSynRatioD\tNonVsSynRatioRegion"); ///////////////////////////////// need modify
//            bw.write("TriadID\tDelVsSynRatioA\tDelVsSynRatioB\tDelVsSynRatioD\tDelVsSynRatioRegion");
            bw.write("TriadID\tDelFreA\tDelFreB\tDelFreD\tDelFreRegion"); //

            bw.newLine();
            int cnt=0;
            int cntremaining =0;
            for (int i = 0; i < tg.getTriadNum(); i++) {
                cnt++;
                String triadID = tg.triadsList.get(i);
                String genea = tg.getGeneinAsub(triadID);
                String geneb = tg.getGeneinBsub(triadID);
                String gened = tg.getGeneinDsub(triadID);

//                String ratioA = genedb.getNonVsSynRatio(genea);
//                String ratioB = genedb.getNonVsSynRatio(geneb);
//                String ratioD = genedb.getNonVsSynRatio(gened);

//                String ratioA = genedb.getDelVsSynRatio(genea);
//                String ratioB = genedb.getDelVsSynRatio(geneb);
//                String ratioD = genedb.getDelVsSynRatio(gened);

                String ratioA = genedb.getPercentageDel(genea);
                String ratioB = genedb.getPercentageDel(geneb);
                String ratioD = genedb.getPercentageDel(gened);


                //filter NA
                if (ratioA.startsWith("N") || ratioB.startsWith("N") || ratioD.startsWith("N")) continue;
                double[] ratiodABD = {Double.parseDouble(ratioA),Double.parseDouble(ratioB),Double.parseDouble(ratioD)};
                String region = Standardization.getNearestPointIndex(ratiodABD).getRegion(); //##引用类来判断所在的位置
                bw.write(triadID+"\t"+ratioA+"\t"+ratioB+"\t"+ratioD+"\t"+region);
                bw.newLine();
                cntremaining++;
            }
            System.out.println(cnt + " triads totally, " + cntremaining + " triads kept");
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    class GeneDB {

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/001_geneSummary_hexaploid.txt.gz";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/002_merge/001_tetraploid_diploid_geneSummary.txt.gz";

        // landrace
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/003_merge/001_LandraceEU_geneSummary.txt.gz";
        //cultivar
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/003_merge/001_Cultivar_geneSummary.txt.gz";

        /**
         * update data
         */
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/002_z1_geneSummary_merge/001_geneSummary_hexaploid.txt.gz";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/102_z1_geneSummary_byChr/001_geneSummary_pseudohexaploid.txt.gz";


        String[] geneArray = AoFile.getgeneArraybyList(infileS,0);
        List<String>[] geneInfo = new List[geneArray.length];

        public GeneDB(){
            this.readFile();
        }

        public String getPercentageHGDeleterious(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(13);
            }else{
                out = "NA";
            }
            return out;
        }

        public String getNonVsSynRatio(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(9);
            }else{
                out = "NA";
            }
            return out;
        }


        public String getDelVsSynRatio(String gene){
            String out = null;
            int synNum = Integer.MIN_VALUE;
            int delNum = Integer.MIN_VALUE;
            double ratio = Double.NaN;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                 synNum = Integer.parseInt(geneInfo[index].get(5));
                 delNum = Integer.parseInt(geneInfo[index].get(12));
                 ratio = (double) delNum/synNum;
                 out = String.valueOf(ratio);
            }else{
                out = "NA";
            }
            return out;
        }

        public String getPercentageDel(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(13);
            }else{
                out = "NA";
            }
            return out;
        }


        public String getPercentageNon(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(8);
            }else{
                out = "NA";
            }
            return out;
        }

        public String getPercentageSyn(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(6);
            }else{
                out = "NA";
            }
            return out;
        }

        public String getIfSiftAligned(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(4);
            }else{
                out = "NA";
            }
            return out;
        }

        public String getSNPNumber(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(2);
            }else{
                out = "NA";
            }
            return out;
        }

        public String getCDSLength(String gene){
            String out = null;
            int index = Arrays.binarySearch(geneArray,gene);
            if (index > -1){
                out = geneInfo[index].get(1);
            }else{
                out = "NA";
            }
            return out;
        }

        public void readFile(){
            AoFile.readheader(infileS);
            Arrays.sort(geneArray);
            for (int i = 0; i < geneArray.length; i++) {
                geneInfo[i] = new ArrayList<>();
            }
            try {
                BufferedReader br = AoFile.readFile(infileS);
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String gene = l.get(0); //去掉转录本
                    gene = gene.split("\\.")[0];
                    int index = Arrays.binarySearch(geneArray,gene);
                    for (int i = 0; i < l.size(); i++) {
                        geneInfo[index].add(l.get(i));
                    }
                }
                br.close();
                System.out.println("Finished read the gene DB");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    /**
     * 将 cultivar 和 landrace_EU 的geneSummary 进行合并
     */
    public void addGroupforgeneSummary(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/003_merge";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/004_addGroup";
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        String group = null;
        String sub = null;
        String outfileS = new File(outfileDirS,"001_LandraceEU_Cultivar_geneSummary.txt.gz").getAbsolutePath();
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(br.readLine() + "\tGroup\tSub");
            bw.newLine();
            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                group = fs[i].getName().split("_")[1];
                br = AoFile.readFile(infileS);
                br.readLine(); //read header
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sub = PStringUtils.fastSplit(temp).get(0).substring(8,9); //第八个字符是基因的亚基因组分类
                    sb.append(temp).append("\t").append(group).append("\t").append(sub);
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

    /**
     * 根据亚群的个体，结合merged VCF 文件，提取亚群的 snp annotation info
     *
     */
    public void getsubspeciesSNPAnnotation(){
//        String taxaList = "";
//        String outfileDirS = "";
        String vcfFileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";

        //需要修改
        //landrace_EU
//        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/009_treeValidatedFroup_byRegion/002_Landrace_European/Landrace_Europe.txt";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/003_hexaploid_subspecies_SNPAnnotation/001_landrace_exon_SNPAnnotation";

        //cultivar
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_bySubspecies/Cultivar.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/003_hexaploid_subspecies_SNPAnnotation/002_cultivar_exon_SNPAnnotation";


        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        AoFile.readheader(fsList.get(0).getAbsolutePath());
        fsList.parallelStream().forEach(f -> {
            String chr = f.getName().substring(0,6); //chr001
            String exonVCF = new File(vcfFileDirS,chr + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            //*********************************** step1: 获取群体的 pos 信息 ***********************************//
            TIntArrayList posl = CalVCF.extractVCFchrPos(exonVCF, taxaList);

            //******* step2: 根据annotation 库， pos列， pols目标集合， 输出文件获取 注释文件的子集 ******//
            String infileS = f.getAbsolutePath();
            //需要修改，输出文件的名字 //需要修改，输出文件的名字 //需要修改，输出文件的名字
//            String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_landraceEU_.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_cultivar_.txt.gz").getAbsolutePath();
            int colmnIndex = 2;
            File out = AoFile.filterTxtLines(infileS,2,posl,outfileS);
        });

    }


    public void mergeTxt(){
//        String infileDirS = "";
//        String outfileS = "";
//        AoFile.mergeTxt(infileDirS,outfileS);

        //合并landrcae_EU的
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/001_Landrcae_EU";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/003_merge/001_LandraceEU_geneSummary.txt";
//        AoFile.mergeTxt(infileDirS,outfileS);


        //合并 cultivar 的
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/002_Cultivar";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/003_merge/001_Cultivar_geneSummary.txt";
//        AoFile.mergeTxt(infileDirS,outfileS);

        //合并 tetraploid diploid
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/001_bychr";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/002_merge/001_tetraploid_diploid_geneSummary.txt.gz";
//        AoFile.mergeTxt(infileDirS,outfileS);

        //合并 hexaploid 的
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/002_geneSummary_byChr";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/002_z1_geneSummary_merge/001_geneSummary_hexaploid.txt.gz";
//        AoFile.mergeTxt(infileDirS,outfileS);

        //合并 pseudohexaploid 的
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/102_geneSummary_byChr";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/102_z1_geneSummary_byChr/001_geneSummary_pseudohexaploid.txt.gz";
        AoFile.mergeTxt(infileDirS,outfileS);


    }

    public void script_getTranscriptSum(){
        //需要修改
//        String infileDirS = ""; // 需要总结的 库 snpAnnotation 文件
//        String outfileDirS =""; // 每条染色体上的基因生成的总结
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/000_hexaploid_SNPAnnotation";
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/001";

        //Landrace_EU
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/003_hexaploid_subspecies_SNPAnnotation/001_landrace_exon_SNPAnnotation";
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/001_Landrcae_EU";

        //Cultivar
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/003_hexaploid_subspecies_SNPAnnotation/002_cultivar_exon_SNPAnnotation";
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/004_geneSummary_byChr/002_Cultivar";

        //Tetraploid
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/001_annotation/001"; // 需要总结的 库 snpAnnotation 文件
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/006_TriadsGeneLoadforPop/002_tetraploid_diploid/002_geneSummary/001_bychr"; // 每条染色体上的基因生成的总结


        /**
         * 2020-11-19 周四 update my data
         */
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/001_exonSNPAnnotation_hexaploid"; // 需要总结的 库 snpAnnotation 文件
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/002_geneSummary_byChr"; // 每条染色体上的基因生成的总结

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/101_exonSNPAnnotation_pseudoHexaploid";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/044_methylation/102_geneSummary_byChr";

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        String[] chrArrAB ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
        String[] chrArrD ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        Arrays.sort(chrArrAB); Arrays.sort(chrArrD);
        for (int j = 0; j < chrArr.length; j++) {

            //#################### 需要修改 ############################//
//            String infileS = "";
//            String outfileS = "";
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_hexaploid_anno.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_hexaploid_geneSummary.txt").getAbsolutePath();

            //Landrace_EU
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_anno_landraceEU_.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_LandraceEU_geneSummary.txt").getAbsolutePath();

            //Cultivar
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_anno_cultivar_.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_Cultivar_geneSummary.txt").getAbsolutePath();

            //Tetraploid and Diploid
//            String infileS = null;
//            String outfileS = null;
//            int index = Arrays.binarySearch(chrArrAB,chrArr[j]);
//            if (index > -1){ // is AABB
//                infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_anno_tetraploid.txt.gz").getAbsolutePath();
//                outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_tetraploid_geneSummary.txt.gz").getAbsolutePath();
//            }else{ //is DD
//                infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_anno_diploid.txt.gz").getAbsolutePath();
//                outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_diploid_geneSummary.txt.gz").getAbsolutePath();
//            }

            // hexaploid
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_hexaploid_anno.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_hexaploid_geneSummary.txt.gz").getAbsolutePath();
//            this.getTranscriptSum_bychr(infileS,outfileS);
//            System.out.println("chr" + chrArr[j] + " is completed at " + outfileS);

            //pseudohexaploid
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_SNP_pseudohexaploid_anno.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_pseudohexaploid_geneSummary.txt.gz").getAbsolutePath();
            this.getTranscriptSum_bychr(infileS,outfileS);
            System.out.println("chr" + chrArr[j] + " is completed at " + outfileS);

        }

    }

    /**
     * 因为多线程运行，posgeneMap超过内存，
     * 故修改思路，进行一条一条计算
     */
    public void getTranscriptSum_bychr(String infileS,String outfileS){
//    public void getTranscriptSum_bychr(){

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/000_hexaploid_SNPAnnotation/chr001_SNP_hexaploid_anno.txt.gz";
        File f = new File(infileS);
        int chrIndex = Integer.parseInt(f.getName().substring(3,6)) -1;
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String outfileS = "/Users/Aoyue/Documents/test.txt";

        double gerpCut = 1;
//        AoFile.readheader(f.getAbsolutePath());
        // 变量命名
        TIntArrayList snpPos = new TIntArrayList();
        TByteArrayList snps = new TByteArrayList();
        TByteArrayList snpAnc = new TByteArrayList();
        List<String> genesList = new ArrayList<>();
        HashMap<Integer, ArrayList<String>> posGeneMap = new HashMap<>(); //1.将posGeneMap建立完整，使每个位点对应哪些基因名字，都装进这个map里
        HashMap<String, Integer> geneCDSLengthMap = new HashMap(); //2.gene <--> cdsLength


        //*********************************** START1 ***********************************//
        /*将所有基因的名字进行for循环输入到数组genes中，对应于每一个基因，我们通过getTranscriptName得到转录本的名字，通过getCDSList方法得到编码序列的长度*/
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chrindex = gf.getGeneChromosome(i)-1;
            if (! (chrindex== chrIndex))continue; //一条一条进行call
            int longTransIndex = gf.getLongestTranscriptIndex(i); //根据基因获取最长转录本索引
            String geneName = gf.getTranscriptName(i, longTransIndex); //根据索引获取最长转录本的名字
            genesList.add(geneName); //将最长转录本的名字加入geneList
            List<Range> cdsList = gf.getCDSList(i, longTransIndex); //根据索引获取最长转录本的CDSList
            int cnt = 0;

            /*对于每一个基因的编码序列，还有很多个cds片段，即cdsList；我们对cdsList进行for循环，得到每个cds的起始和终止位置，从而计算出总长*/
            for (int j = 0; j < cdsList.size(); j++) {
                int rStart = cdsList.get(j).start;
                int rEnd = cdsList.get(j).end;
                for (int k = rStart; k < rEnd; k++) {
                    /*posGeneMap是一个HashMap数组，一条染色体对应一个String类型的ArrayList；
                    故得到该位点的所属基因名字列表，如果该位点不含基因名，就将 genename赋值给该位点，完善posGeneMap
                    否则，如果该位点含有其他基因的基因名字，依旧把genename赋值给该位点*/
                    ArrayList<String> geneNameList = posGeneMap.get(k); //建立map的关系，那个位点对应哪个list HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
                    if (geneNameList == null) {
                        geneNameList = new ArrayList();
                        geneNameList.add(geneName);
                        posGeneMap.put(k, geneNameList);
                    }
                    else {
                        geneNameList.add(geneName);
                        posGeneMap.put(k, geneNameList); /*最终将posGeneMap绘图完成*/
                    }
                    cnt++; /*每一个CDS位点相加，最终得到这个cds的长度。*/
                } //k end
                // 最终cnt是一个基因的所有cdslist中，每个cds的每个位点包含的基因数目的总和
            } //该循环是一个基因的所有cds循环
            geneCDSLengthMap.put(geneName, cnt); //
        } //gene end

        System.out.println(genesList.size()  + " genes in the pgf file");
        //*********************************** END1 ***********************************//

        //*********************************** START 2 ***********************************//
        //
        String[] genes = genesList.toArray(new String[genesList.size()]);
        Arrays.sort(genes); //genes指所有基因对应最长转录本的名字的组合，是一个数组。
        int[] snpCount = new int[genes.length];

        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length]; //非同义突变的数目
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] nasiftCount = new int[genes.length];  //非同义突变，但是sift值是NA的数目
        TIntArrayList delPosList = new TIntArrayList(); //有害变异的位点集合

        int[] CSSynCount = new int[genes.length];
        int[] CSNonCount = new int[genes.length];
        int[] CSDelCount = new int[genes.length];
        int[] CSDelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];

        int[] gerpAlignCount = new int[genes.length];
        double[] gerpScore = new double[genes.length];

            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            TByteArrayList snpAncList = new TByteArrayList();
            try{
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++; //检测程序运行的情况
                    if (cnt%1000000 == 0) System.out.println("chr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ###hmpInfo Process");
                    l = PStringUtils.fastSplit(temp);

//                    int pos = Integer.valueOf(l.get(2));
//                    String ref = l.get(3);
//                    String alt = l.get(4);
//                    String ancestral = l.get(31);
//                    String variant_type = l.get(12);
//                    String sift = l.get(13);
//                    String trans = l.get(10);
//                    String gerp = l.get(20);

                    int pos = Integer.valueOf(l.get(2));
                    String ref = l.get(3);
                    String alt = l.get(4);
                    String ancestral = l.get(15);
                    String variant_type = l.get(12);
                    String sift = l.get(16);
                    String trans = l.get(10);
                    String gerp = l.get(20);

                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap.get(pos); //根据pos信息，得到该pos对应的gene name的集合
                    if (geneNameList == null) continue; //说明该变异位点不在基因区
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        int geneIndex = Arrays.binarySearch(genes, gene);//在基因库的第i个位置，该基因数加一
                        snpCount[geneIndex]++; //第i个位置的基因含有的snp数目；

                        if (gerp.startsWith("NA"))continue;
                        double scoreValue = Double.parseDouble(gerp);
                        gerpAlignCount[geneIndex]++; //如果枝长不是0，说明该位点存在保守不保守
                        gerpScore[geneIndex]+=scoreValue; //第i个基因的gerpscore的总和是多少
                    }
                    /*将该位置加入snpPosList，表明该位置是有基因的，没有基因的位置就不必加到snpPosList中去*/
                    snpPosList.add(pos);
                    snpList.add(alt.getBytes()[0]); //返回该位点的alt碱基的二进制码
                    snpAncList.add(ancestral.getBytes()[0]);  //返回该位点的ancestral碱基的二进制码

                    if (variant_type.equals("NA")) continue;
                    if (trans.equals("NA")) continue;
                    int geneIndex = Arrays.binarySearch(genes, trans); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) continue;
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (ancestral.equals(alt)) { //如果ancestral allele 等于alt的话，derived allele就等于1
                        derivedState = 1; //mean cs carries derived allele
                    }
                    else if (ancestral.equals(ref)) { //如果ancestral allele 等于ref的话，derived allele就等于0
                        derivedState = 0;
                    }

                    if (derivedState == -1) noAncCount[geneIndex]++; //如果ancestral allele 不存在的话，derived allele就等于-1

                    String  type = null;
                    if (variant_type.equals("NA"))continue;
                    if (variant_type.equals("SYNONYMOUS")){
                        type = "Syn";
                        synCount[geneIndex]++;
                        if (derivedState == 1){
                            CSSynCount[geneIndex]++;
                        }
                    }
                    if (variant_type.equals("NONSYNONYMOUS")){
                        type = "Non";
                        nonCount[geneIndex]++;
                        if (derivedState == 1){
                            CSNonCount[geneIndex]++;
                        }

                        if (sift.startsWith("NA")){
                            nasiftCount[geneIndex]++; //是非同义突变，但是没有sift值
                        }
                        else if (!sift.startsWith("NA")){
                            double siftd = Double.parseDouble(sift);
                            if (siftd < 0.05){ //start1
                                delCount[geneIndex]++;
                                delPosList.add(pos);
                                if (derivedState == 1) { //说明 derived allele是 ref
                                    CSDelCount[geneIndex]++;
                                }

                                if (gerp.startsWith("NA"))continue;
                                double scoreValue = Double.parseDouble(gerp);
                                if (scoreValue >= gerpCut) {
                                    delHGCount[geneIndex]++;
                                }
                                if (derivedState == 1) {
                                    CSDelHGCount[geneIndex]++;
                                }
                            } //end1

                        }
                    } //nonsynonymous 结束
                } //while end
                br.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        //*********************************** END2 ***********************************//

        //*********************************** START 3 ***********************************//
        //

        gf.sortGeneByStartPosition();
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpScore";
            header = header +"\tNumAmbigousAnc\tCSNumberOfSyn\tCSPercentageSyn\tCSNumberOfNon\tCSPercentageNon\tCSNumberOfDeleterious\tCSPercentageDeleterious\tCSNumberOfHGDeleterious\tCSPercentageHGDeleterious";
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < genes.length; i++) {
                StringBuilder sb =new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append(String.format("%.4f",(double)snpCount[i]/cdsLength)).append("\t");
                int ifSiftAligned = 1;
                if (nasiftCount[i] == nonCount[i]) ifSiftAligned = 0; //非同义突变，但是sift值是NA的数目 等于 非同义突变的数目，那么说明在这个基因内部没有sift突变
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append(String.format("%.4f",(double)synCount[i]/cdsLength)).append("\t");
                sb.append(nonCount[i]).append("\t").append(String.format("%.4f",(double)nonCount[i]/cdsLength)).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) ratio = Double.NaN;
                else ratio = (double)nonCount[i]/synCount[i];
                sb.append(String.format("%.4f",ratio)).append("\t").append(delCount[i]).append("\t").append(String.format("%.4f",(double)delCount[i]/cdsLength)).append("\t");
                sb.append(delHGCount[i]).append("\t").append(String.format("%.4f",(double)delHGCount[i]/cdsLength)).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) ifGerpAligned = 0;
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append(String.format("%.4f",(double)gerpAlignCount[i]/cdsLength)).append("\t");

                if (gerpAlignCount[i] == 0) sb.append(Double.NaN);
                else sb.append(String.format("%.4f",(double)gerpScore[i]/gerpAlignCount[i]));

//                if (snpGerpAlignCount[i] == 0) sb.append(Double.NaN);
//                else sb.append((double)snpGerpScore[i]/snpGerpAlignCount[i]);

                double cdsL = (double)(snpCount[i]-noAncCount[i])/snpCount[i]*cdsLength; //有祖先状态的snp除以总共snp乘以cds长度
                sb.append("\t").append(noAncCount[i]).append("\t").append(CSSynCount[i]).append("\t").append(String.format("%.4f",(double)CSSynCount[i]/cdsL)).append("\t");
                sb.append(CSNonCount[i]).append("\t").append(String.format("%.4f",(double)CSNonCount[i]/cdsL)).append("\t");

                sb.append(CSDelCount[i]).append("\t").append(String.format("%.4f",(double)CSDelCount[i]/cdsL)).append("\t");
                sb.append(CSDelHGCount[i]).append("\t").append(String.format("%.4f",(double)CSDelHGCount[i]/cdsL));

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

    /**
     * 获取转录本的总结
     */
    public void getTranscriptSum(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/001/transcriptSummary.txt";

        double gerpCut = 1;
        List<File> fs = AoFile.getFileListInDir(infileDirS);
        AoFile.readheader(fs.get(0).getAbsolutePath());
        int chrNum = fs.size();
        // 变量命名
        int[][] snpPos = new int[chrNum][]; //每条染色体的snp pos位点
        byte[][] snps = new byte[chrNum][]; //每条染色体pos对应的ALT碱基
        byte[][] snpAnc = new byte[chrNum][]; //每条染色体pos对应的 ancestral 碱基

        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        List<String> genesList = new ArrayList<>();
        HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum]; //1.将posGeneMap建立完整，使每个位点对应哪些基因名字，都装进这个map里
        for (int i = 0; i < chrNum; i++) {
            posGeneMap[i] = new HashMap();
        }
        HashMap<String, Integer> geneCDSLengthMap = new HashMap(); //2.gene <--> cdsLength

        //*********************************** START1 ***********************************//
        /*将所有基因的名字进行for循环输入到数组genes中，对应于每一个基因，我们通过getTranscriptName得到转录本的名字，通过getCDSList方法得到编码序列的长度*/
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chrIndex = gf.getGeneChromosome(i)-1;
            if ((chrIndex== -1))continue;
            int longTransIndex = gf.getLongestTranscriptIndex(i); //根据基因获取最长转录本索引
            String geneName = gf.getTranscriptName(i, longTransIndex); //根据索引获取最长转录本的名字
            genesList.add(geneName); //将最长转录本的名字加入geneList
            List<Range> cdsList = gf.getCDSList(i, longTransIndex); //根据索引获取最长转录本的CDSList
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
                    }
                    else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList); /*最终将posGeneMap绘图完成*/
                    }
                    cnt++; /*每一个CDS位点相加，最终得到这个cds的长度。*/
                }

                // 最终cnt是一个基因的所有cdslist中，每个cds的每个位点包含的基因数目的总和
            } //该循环是一个基因的所有cds循环
            geneCDSLengthMap.put(geneName, cnt); //
        }

        System.out.println(genesList.size()  + " genes in the pgf file");
        //*********************************** END1 ***********************************//

        //*********************************** START 2 ***********************************//
        //
        String[] genes = genesList.toArray(new String[genesList.size()]);
        Arrays.sort(genes); //genes指所有基因对应最长转录本的名字的组合，是一个数组。
        int[] snpCount = new int[genes.length];

        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length]; //非同义突变的数目
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];  //非同义突变，但是sift值是NA的数目
        TIntArrayList[] delPosList = new TIntArrayList[chrNum]; //有害变异的位点集合

        int[] CSSynCount = new int[genes.length];
        int[] CSNonCount = new int[genes.length];
        int[] CSDelCount = new int[genes.length];
        int[] CSDelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];

        int[] gerpAlignCount = new int[genes.length];
        int[] snpGerpAlignCount = new int[genes.length];
        double[] gerpScore = new double[genes.length];
        double[] snpGerpScore = new double[genes.length];



        fs.stream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3,6)) -1;
            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            TByteArrayList snpAncList = new TByteArrayList();
            delPosList[chrIndex] = new TIntArrayList();
            try{
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++; //检测程序运行的情况
                    if (cnt%1000000 == 0) System.out.println("chr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ###hmpInfo Process");
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(2));
                    String ref = l.get(3);
                    String alt = l.get(4);
                    String ancestral = l.get(31);
                    String variant_type = l.get(12);
                    String sift = l.get(13);
                    String trans = l.get(10);
                    String gerp = l.get(20);



                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos); //根据pos信息，得到该pos对应的gene name的集合
                    if (geneNameList == null) continue; //说明该变异位点不在基因区
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        int geneIndex = Arrays.binarySearch(genes, gene);//在基因库的第i个位置，该基因数加一
                        snpCount[geneIndex]++; //第i个位置的基因含有的snp数目；

                        if (gerp.startsWith("NA"))continue;
                        double scoreValue = Double.parseDouble(gerp);
                        gerpAlignCount[geneIndex]++; //如果枝长不是0，说明该位点存在保守不保守
                        gerpScore[geneIndex]+=scoreValue; //第i个基因的gerpscore的总和是多少
                        snpGerpAlignCount[geneIndex]++;
                        snpGerpScore[geneIndex]+=scoreValue;
                    }
                    /*将该位置加入snpPosList，表明该位置是有基因的，没有基因的位置就不必加到snpPosList中去*/
                    snpPosList.add(pos);
                    snpList.add(alt.getBytes()[0]); //返回该位点的alt碱基的二进制码
                    snpAncList.add(ancestral.getBytes()[0]);  //返回该位点的ancestral碱基的二进制码


                    if (variant_type.equals("NA")) continue;
                    if (trans.equals("NA")) continue;
                    int geneIndex = Arrays.binarySearch(genes, trans); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) continue;
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (ancestral.equals(alt)) { //如果ancestral allele 等于alt的话，derived allele就等于1
                        derivedState = 1; //mean cs carries derived allele
                    }
                    else if (ancestral.equals(ref)) { //如果ancestral allele 等于ref的话，derived allele就等于0
                        derivedState = 0;
                    }

                    if (derivedState == -1) noAncCount[geneIndex]++; //如果ancestral allele 不存在的话，derived allele就等于-1

                    String  type = null;
                    if (variant_type.equals("NA"))continue;
                    if (variant_type.equals("SYNONYMOUS")){
                        type = "Syn";
                        synCount[geneIndex]++;
                        if (derivedState == 1){
                            CSSynCount[geneIndex]++;
                        }
                    }
                    if (variant_type.equals("NONSYNONYMOUS")){
                        type = "Non";
                        nonCount[geneIndex]++;
                        if (derivedState == 1){
                            CSNonCount[geneIndex]++;
                        }
                        if (sift.startsWith("NA")){
                            naCount[geneIndex]++; //是非同义突变，但是没有sift值
                        }
                        else if (!sift.startsWith("NA")){
                            double siftd = Double.parseDouble(sift);
                            if (siftd < 0.05){ //start1
                                delCount[geneIndex]++;
                                delPosList[chrIndex].add(pos);
                                if (derivedState == 1) { //说明 derived allele是 ref
                                    CSDelCount[geneIndex]++;
                                }

                                if (gerp.startsWith("NA"))continue;
                                double scoreValue = Double.parseDouble(gerp);
                                if (scoreValue > gerpCut) {
                                    delHGCount[geneIndex]++;
                                }
                                if (derivedState == 1) {
                                    CSDelHGCount[geneIndex]++;
                                }
                            } //end1

                        }
                    } //nonsynonymous 结束
                } //while end
                br.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            snpPos[chrIndex] = snpPosList.toArray(); //第n条染色体的在基因区间的snp所对应的pos的集合；
            snps[chrIndex] = snpList.toArray(); //第n条染色体的在区间内的snp所对应的alt的集合
            snpAnc[chrIndex] = snpAncList.toArray(); //第n条染色体的在区间内的snp所对应的ancestral allele的集合

//            int[][] delPos = new int[chrNum][];
//            for (int i = 0; i < chrNum; i++) {
//                delPos[i] = delPosList[i].toArray();
//                Arrays.sort(delPos[i]);
//            }
        });

        //*********************************** END2 ***********************************//

        //*********************************** START 3 ***********************************//
        //

        gf.sortGeneByStartPosition();
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpScore\tSNPMeanGerpScore";
            header = header +"\tNumAmbigousAnc\tCSNumberOfSyn\tCSPercentageSyn\tCSNumberOfNon\tCSPercentageNon\tCSNumberOfDeleterious\tCSPercentageDeleterious\tCSNumberOfHGDeleterious\tCSPercentageHGDeleterious";
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < genes.length; i++) {
                StringBuilder sb =new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append((double)snpCount[i]/cdsLength).append("\t");
                int ifSiftAligned = 1;
                if (naCount[i] == nonCount[i]) ifSiftAligned = 0; //非同义突变，但是sift值是NA的数目 等于 非同义突变的数目，那么说明在这个基因内部没有sift突变
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append((double)synCount[i]/cdsLength).append("\t");
                sb.append(nonCount[i]).append("\t").append((double)nonCount[i]/cdsLength).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) ratio = Double.NaN;
                else ratio = (double)nonCount[i]/synCount[i];
                sb.append(ratio).append("\t").append(delCount[i]).append("\t").append((double)delCount[i]/cdsLength).append("\t");
                sb.append(delHGCount[i]).append("\t").append((double)delHGCount[i]/cdsLength).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) ifGerpAligned = 0;
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append((double)gerpAlignCount[i]/cdsLength).append("\t");

                if (gerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t");
                else sb.append((double)gerpScore[i]/gerpAlignCount[i]).append("\t");

                if (snpGerpAlignCount[i] == 0) sb.append(Double.NaN);
                else sb.append((double)snpGerpScore[i]/snpGerpAlignCount[i]);

                double cdsL = (double)(snpCount[i]-noAncCount[i])/snpCount[i]*cdsLength;
                sb.append("\t").append(noAncCount[i]).append("\t").append(CSSynCount[i]).append("\t").append((double)CSSynCount[i]/cdsL).append("\t");
                sb.append(CSNonCount[i]).append("\t").append((double)CSNonCount[i]/cdsL).append("\t");

                sb.append(CSDelCount[i]).append("\t").append((double)CSDelCount[i]/cdsL).append("\t");
                sb.append(CSDelHGCount[i]).append("\t").append((double)CSDelHGCount[i]/cdsL);

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

    public void getHexaploidAnnotation(){
        String inputDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";
        String outDir = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/000_hexaploid_SNPAnnotation";
        List<File> files= AoFile.getFileListInDir(inputDir);
        String[] outNames=files.stream().map(File::getName).map(s->s.replaceAll("_anno.txt.gz", "_hexaploid_anno.txt")).toArray(String[]::new);
        RowTableTool<String> rowTable;
        Predicate<List<String>> removed= l->Double.parseDouble(l.get(8))==0;
        for (int i = 0; i < files.size(); i++) {
            rowTable=new RowTableTool<>(files.get(i).getAbsolutePath());
            rowTable.removeIf(removed);
            rowTable.write(new File(outDir, outNames[i]).getAbsolutePath());

        }
    }
}
