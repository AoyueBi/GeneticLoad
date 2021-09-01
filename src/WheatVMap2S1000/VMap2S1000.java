package WheatVMap2S1000;

import AoUtils.*;
import AoUtils.Gene.GeneMisc;
import PopulationAnalysis.XPCLR;
import WheatGeneticLoad.FilterVCF2;
import WheatGeneticLoad.VariantsSum;
import daxing.common.IOTool;
import daxing.common.PGF;
import daxing.common.RowTableTool;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaByte;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import pgl.infra.window.SimpleWindow;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class VMap2S1000 {
    public VMap2S1000(){
        /**
         * Finalize VMap 2.0_S1062
         */
//        this.bgzip();
//        this.rename();
//        this.filterSNPtoBi_parallel();
//        this.filterN_fromVCF();


        /**
         * vcf QC
         */
//        new FilterVCF2().QC();
//        new FilterVCF2().statVcfDepth_SD();
//        new FilterVCF2().mkDepthOfVMapII();
//        new FilterVCF2().mkDepthSummary();
//        this.extractVCF(); //提取不同倍性的VCF
//        new FilterVCF2().QC();

//        this.getVCFbyPolyploid();
//        this.getMAFcountfromPop();

        /**
         * gene site annotation
         */
//        this.snpAnnotationBuild(); //include many methods XXXXXXX
//        new DeleteriousCount();

        /**
         * XPCLR
         */
//        new XPCLR();
//        new DeleteriousXPCLRS1000();
//        this.getVCF(); //提取亚群的VCF文件
//        new GeneMisc();

        /**
         *  Rht gene 的验证
         */
//        this.getVCFofRht("Rht-B1");
//        this.getVCFofRht("Rht-D1");
//        new Model().runJarParallele(); //checkAltCaseCount()
//        this.getGenoTable();

        /**
         * 全基因组 del 分布
         */
//        new VariantsSum().AddGenePosition();
        this.WindowDelvsSyn_fromExonAnnotation();


        /**
         * 有害突变数目随样本大小变化
         */

    }


    /**
     * 这里只研究 del 和 syn
     */
    public void WindowDelvsSyn_fromExonAnnotation(){
//        int windowSize = 2000000; //2 M
//        int windowStep = 1000000; //1 M

//        int windowSize = 20000000; //20 M
//        int windowStep = 5000000; //5 M

        int windowSize = 10000000; //2 M
        int windowStep = 1000000; //1 M


//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_" + windowSize + "Window" + windowStep + "step.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_" + windowSize + "Window" + windowStep + "step.txt";

        //******* 2021-09-01 周三 *******//
        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/001_geneSNPAnno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/009_genomeScan_delvcSyn/001/001_delVSsynOnChr_" + windowSize + "window" + windowStep + "step.txt";
        String outfile2S = outfileS.replaceFirst(".txt","_addEffectiveCDSLength.txt");
        AoFile.readheader(infileS);
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int chrNum = chrArr.length;

        int[][] synPos = new int[chrNum][];
        int[][] nonsynPos = new int[chrNum][];
        int[][] siftGerpPos = new int[chrNum][];
        int[][] siftPos = new int[chrNum][];
        int[][] gerpPos = new int[chrNum][];
        int[][] stopGainSiftPos = new int[chrNum][];
        int[][] vepPos = new int[chrNum][];
        int[][] stopGainVEPPos = new int[chrNum][];
        int[][] snpEffPos = new int[chrNum][];

        TIntArrayList[] synposList = new TIntArrayList[chrNum];
        TIntArrayList[] nonsynPosList = new TIntArrayList[chrNum];
        TIntArrayList[] siftGerpPosList = new TIntArrayList[chrNum];
        TIntArrayList[] siftPosList = new TIntArrayList[chrNum];
        TIntArrayList[] gerpPosList = new TIntArrayList[chrNum];
        TIntArrayList[] stopGainSiftPosList = new TIntArrayList[chrNum];
        TIntArrayList[] vepPosList = new TIntArrayList[chrNum];
        TIntArrayList[] stopGainVEPPosList = new TIntArrayList[chrNum];
        TIntArrayList[] snpEffPosList = new TIntArrayList[chrNum];


        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            synposList[i] = new TIntArrayList();
            nonsynPosList[i] = new TIntArrayList();
            siftGerpPosList[i] = new TIntArrayList();
            siftPosList[i] = new TIntArrayList();
            gerpPosList[i] = new TIntArrayList();
            stopGainSiftPosList[i] = new TIntArrayList();
            vepPosList[i] = new TIntArrayList();
            stopGainVEPPosList[i] = new TIntArrayList();
            snpEffPosList[i] = new TIntArrayList();
        }

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr =  Integer.parseInt(l.get(1));
                int pos = Integer.parseInt(l.get(2));

                String chromosome = RefV1Utils.getChromosome(chr,1); //获取 RefChr
                int posonchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                int index = Arrays.binarySearch(chrArr,chromosome); //染色体的索引号
                String ancestralAllele = l.get(9);
                String variantType = l.get(13);
                String sift = l.get(16); // derived_sift
                String gerp = l.get(11);
                String Effect_VEP = l.get(17);
                String Impact_VEP = l.get(18);
                String Impact_snpEff = l.get(20);

                // 只考虑有 ancestral 状态的那些位点
                String ref = l.get(3); String alt = l.get(4);
                if (!ancestralAllele.equals(ref) && !ancestralAllele.equals(alt)) continue;

                if (variantType.equals("SYNONYMOUS")){
                    synposList[index].add(posonchromosome);
                }

                if (variantType.equals("NONSYNONYMOUS")){
                    nonsynPosList[index].add(posonchromosome); // nonsyn
                }

                if (variantType.equals("NONSYNONYMOUS")){ // sift_gerp
                    if (!sift.startsWith("N")){
                        if(!gerp.startsWith("N")){
                            double gerpd = Double.parseDouble(gerp);
                            double siftd = Double.parseDouble(sift);
                            if (siftd < 0.05 && gerpd > 1){
                                siftGerpPosList[index].add(posonchromosome);
                            }
                        }
                    }
                }

                if (variantType.equals("NONSYNONYMOUS")){ // sift
                    if (!sift.startsWith("N")){
                        double siftd = Double.parseDouble(sift);
                        if (siftd < 0.05){
                            siftPosList[index].add(posonchromosome);
                        }
                    }
                }

                if (variantType.equals("NONSYNONYMOUS")){ // gerp
                    if(!gerp.startsWith("N")){
                        double gerpd = Double.parseDouble(gerp);
                        if (gerpd > 1){
                            gerpPosList[index].add(posonchromosome);
                        }
                    }
                }

                if (variantType.equals("STOP-GAIN")){ // stop gain_sift
                    stopGainSiftPosList[index].add(posonchromosome);
                }

                if (Impact_VEP.equals("HIGH")){ // vepPosList
                    vepPosList[index].add(posonchromosome);
                }

                if (Effect_VEP.contains("start_lost") || Effect_VEP.contains("stop_gained") || Effect_VEP.contains("stop_lost")){ // vepPosList
                    stopGainVEPPosList[index].add(posonchromosome);
                }

                if (Impact_snpEff.equals("HIGH")){ // vepPosList
                    snpEffPosList[index].add(posonchromosome);
                }
            }
            System.out.println("======== completing the posList DB on all chromosome.");

            String[] variantTypeArray = {"001_synonymous","002_nonsynonymous","003_nonsynGERPandDerivedSIFT","004_nonsynDerivedSIFT","005_GERP",
            "006_StopGain","007_VEP","008_snpEff","009_VEP_stopGained"};
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String outheader = "CHROM\tBIN_START\tBIN_END\tBIN_START_scale"; //和vcftools的格式保持一致
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < variantTypeArray.length ; i++) {
                sb.append("\t").append(variantTypeArray[i]).append("_Count");
            }
            bw.write(outheader);bw.write(sb.toString());bw.newLine();

            sb.setLength(0);
            for (int i = 0; i < chrArr.length; i++) {
                String chromosome = chrArr[i];
                int chrLength = RefV1Utils.getChromosomeLength(chromosome);
                SimpleWindow sw = new SimpleWindow(chrLength, windowSize, windowStep);

                sw.addPositionCount(synposList[i].toArray());
                int[] synPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(nonsynPosList[i].toArray());
                int[] nonsynPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(siftGerpPosList[i].toArray());
                int[] siftGerpPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(siftPosList[i].toArray());
                int[] siftPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(gerpPosList[i].toArray());
                int[] gerpPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(stopGainSiftPosList[i].toArray());
                int[] stopGainSiftPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(vepPosList[i].toArray());
                int[] vepPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(stopGainVEPPosList[i].toArray());
                int[] stopGainVEPPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(snpEffPosList[i].toArray());
                int[] snpEffPosWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                int[] windowStarts = sw.getWindowStarts();
                int[] windowEnds = sw.getWindowEnds();
                for (int j = 0; j < windowStarts.length; j++) {
                    sb.setLength(0);
                    String posscale = WheatUtils.getScaledPos(chromosome,windowStarts[j]);
                    sb.append(chromosome).append("\t").append(windowStarts[j]).append("\t").append(windowEnds[j]).append("\t").append(posscale).append("\t");
                    sb.append(synPosWindowCount[j]).append("\t").append(nonsynPosWindowCount[j]).append("\t");
                    sb.append(siftGerpPosWindowCount[j]).append("\t").append(siftPosWindowCount[j]).append("\t");
                    sb.append(gerpPosWindowCount[j]).append("\t").append(stopGainSiftPosWindowCount[j]).append("\t");
                    sb.append(vepPosWindowCount[j]).append("\t").append(stopGainVEPPosWindowCount[j]).append("\t");
                    sb.append(snpEffPosWindowCount[j]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println("======== " + chromosome + " is completed.");
            }
            bw.flush();bw.close();
//*********************************** add cds length *******************************************
            this.addCDSLengthInWindow(outfileS, windowSize, windowStep,outfile2S);
            System.out.println("======== " + "add effective CDS length.");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 将 del nonsyn syn 的结果写在同一个文件中，不进行分开， 和addCDSLengthInWindow方法的区别在于，前者只计算了del 和 syn 的count 和cds
     * @param dbFileS
     * @param windowSize
     * @param windowStep
     * @param outfile2S
     */
    public void addCDSLengthInWindow (String dbFileS, int windowSize, int windowStep, String outfile2S) {
//        int windowSize = 2000000;
//        int windowStep = 1000000;

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String dbFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_2000000Window1000000step.txt";
//        String hcGeneFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        String nonoverlapGeneFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";

        RowTable<String> gt = new RowTable<>(nonoverlapGeneFileS);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName(); //通过名字排序
        List<String> chromosomeList = RefV1Utils.getChromosomeList(); //返回1A 2A - 7D
        TIntArrayList cdsWindowList = new TIntArrayList(); //所有cds的list，这里并没有指明大小

        for (int i = 0; i < chromosomeList.size(); i++) { //每条染色体进行循环
            // 研究思路： 首先得到每条染色体的长度，建立一个 simpleWindow 类，
            // 对基因列表中的所有基因进行循环，如果输入的基因染色体不在该循环染色体内，就跳出循环。
            //根据基因名字，获取基因的Index，后续找到最长转录本，获取该转录本的 cdsList ,对 cdsList 进行循环，一一装入 SimpleWindow 类中，调用 addPositionCountFromRange 方法，统计不同 Bin 里的 cds 长度。
            int chrlength = RefV1Utils.getChromosomeLength(chromosomeList.get(i));
            SimpleWindow sw = new SimpleWindow(chrlength, windowSize, windowStep); //new 一个 SimpleWindow 类
            int chrID = RefV1Utils.getChrID(chromosomeList.get(i), 1); //根据chromosomeList中的1A等，获取chrID
            int geneIndex = -1;
            int tranIndex = -1;
            int cdsStart = -1;
            int cdsEnd = -1;
            for (int j = 0; j < gt.getRowNumber(); j++) { //gt已通过名字进行了排序
                int isuniqueGene = Integer.parseInt(gt.getCell(j, 3));
                if (isuniqueGene == 0) continue;
                int currentChrID = Integer.parseInt(gt.getCell(j, 5));
                if (currentChrID < chrID) continue; //如果，gt文件中的chrID小于循环内的chromosomeList.get(i)，则跳出循环，一直到gt列表里的chrID等于此时正在循环的染色体号
                else if (currentChrID > (chrID +1)) break; //如果，gt文件中的currentChrID 等于chrID+1，说明还是在1A内，其他情况都终止循环
                geneIndex = gf.getGeneIndex(gt.getCell(j, 0)); //根据基因名字（不是转录本的名字，没有.后缀12）获取索引
                for (int k = 0; k < gf.getTranscriptNumber(geneIndex); k++) { //获取最长转录本的index
                    if (!gt.getCell(j,4).equals(gf.getTranscriptName(geneIndex, k)))continue; //如果 gt.getCell(j,1) 是该基因的最长转录本，已经提前总结出
                    tranIndex = k; //该循环的意思是：如果最长转录本不等于k，那么就跳出循环，一直到等于k为止，最后终止循环。
                    break;
                }
                List<Range> cdsList = gf.getCDSList(geneIndex, tranIndex); //根据 gene index 和 tranIndex 获取该转录本的 cdsList
                for (int k = 0; k < cdsList.size(); k++) {
                    cdsStart = cdsList.get(k).start;
                    cdsEnd = cdsList.get(k).end;
                    if (currentChrID == chrID + 1) { //chrID 是根据chromosomeList求出来的1A的第一个染色体chr001,同时chr002也属于1A。 这里如果是在2号的话，则需要改变位置。
                        cdsStart = RefV1Utils.getPosOnChromosome(currentChrID, cdsStart);
                        cdsEnd = RefV1Utils.getPosOnChromosome(currentChrID, cdsEnd);
                    }
                    sw.addPositionCountFromRange(cdsStart, cdsEnd); //每个基因的每个cds的起始和终止段，都加入了 sw类中
                }
            } //所有的基因都进行了cds的起始和终止的计数，并且加入了 sw类中
            cdsWindowList.add(sw.getWindowValuesInt()); //将每个window的值装入cdsWindowList中
        }


        try {
            BufferedReader br = AoFile.readFile(dbFileS);
            //dbFileS include: CHROM	BIN_START	BIN_END	BIN_START_scale	DelCount	NonsynCount	SynCount	DelSynRatio	NonsynSynRatio
            BufferedWriter bw = AoFile.writeFile(outfile2S);
            String header = br.readLine();
            bw.write(header + "\tCDSLength");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int line = 0;
            int cdsLength = -1;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                cdsLength = cdsWindowList.get(line);
                line++;
                sb.append(temp).append("\t").append(cdsLength);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println( "======== completed at " + outfile2S);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void getGenoTable(){
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF/chr021_RhtB1_vmap2.1.vcf.gz";
//        String taxaListFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/taxaList/LR_EU_CL.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/chr021_RhtB1_genoTable.txt.gz";

//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF/chr023_RhtD1_vmap2.1.vcf.gz";
//        String taxaListFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/taxaList/LR_EU_CL.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/chr023_RhtD1_genoTable.txt.gz";


//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF/chr021_RhtB1_vmap2.1.vcf.gz";
//        String taxaListFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/taxaList/TaxaList_HexaploidandTetraploid.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/S1026/chr021_RhtB1_genoTable.txt.gz";

        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF/chr023_RhtD1_vmap2.1.vcf.gz";
        String taxaListFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/taxaList/TaxaList_HexaploidandDiploid.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/002_genoTable/S1026/chr023_RhtD1_genoTable.txt.gz";


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
//                    bw.write("CHROM" + "\t" + l.get(1));
                    bw.write("CHROM" + "\t" + l.get(1) + "\t" + l.get(3) + "\t" + l.get(4));
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
//                    int nAlt = PStringUtils.fastSplit(altList, ",").size();
//                    if (nAlt > 1) continue; //filter alt num with 2 or more
//                    if (altList.equals("D") || altList.equals("I")) continue; //filter D I

                    String chr = l.get(0);
                    String pos = l.get(1);
                    String ref = l.get(3);
                    String alt = l.get(4);
//                    bw.write(chr + "\t" + pos);
                    bw.write(chr + "\t" + pos + "\t" + ref + "\t" + alt);

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
     *
     */
    public void checkAltCaseCount(String infileS, String outfileS){
//        String infileS = "/Users/Aoyue/Documents/chr021_RhtB1_vmap2.1.vcf.gz";
//        String outfileS ="/Users/Aoyue/Documents/countCase.txt";
        Set<String> altSet = new HashSet<>();
        List<String> altList = new ArrayList<>();
        int cnt = 0;
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                if (temp.startsWith("#"))continue;
                l = PStringUtils.fastSplit(temp);
                String alt = l.get(4);
                altSet.add(alt);
                altList.add(alt);
                cnt++;
            }

            br.close();bw.flush();bw.close();



        }catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }


        File out = new File(outfileS);
        BufferedWriter bw = AoFile.writeFile(out.getAbsolutePath());
        System.out.println(altList.size() + " list个数");
        Set<String> s = new HashSet<>(altList);
        System.out.println(s.size() + " set个数");
        System.out.println(s);

        try{
            bw.write("Case\tCount");
            bw.newLine();
            for(String a : s){
                System.out.println(a + "\t" + Collections.frequency(altList, a));
                bw.write(a+"\t"+Collections.frequency(altList, a));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void getVCFofRht(String geneGoal){
        String infileS = null;
        String outfileS =null;
        String[] geneArray = {"Rht-B1","Rht-D1"};


        String gene = null;
        gene = geneGoal;

        // // RhtB1 position
        int chrID4b = 21;
        int pos4b_1 = 30861382;
        int pos4b_2 = 30863283;
        // // RhtD1 position
        int chrID4d = 23;
        int pos4d_1 = 18781062;
        int pos4d_2 = 18782934;

//        String infileS1 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/007_geneVCF/chr021_gene_vmap2.1.vcf.gz";
//        String outfileS1 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF/chr021_Rht_vmap2.1.vcf.gz";
//        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/007_geneVCF/chr023_gene_vmap2.1.vcf.gz";
//        String outfileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF/chr023_Rht_vmap2.1.vcf.gz";

        String infileS1 = "/data1/home/xinyue/Vmap2_Out/chr021.vmap2.vcf.gz";
        String outfileS1 = "/data4/home/aoyue/vmap2/analysis/042_Rht/001_outVCF/chr021_RhtB1_vmap2.1.vcf.gz";
        String infileS2 = "/data1/home/xinyue/Vmap2_Out/chr023.vmap2.vcf.gz";
        String outfileS2 = "/data4/home/aoyue/vmap2/analysis/042_Rht/001_outVCF/chr023_RhtD1_vmap2.1.vcf.gz";

        int chrRht = Integer.MIN_VALUE;
        int pos1Rht = Integer.MIN_VALUE;
        int pos2Rht = Integer.MIN_VALUE;
        if (gene.equals(geneArray[0])){
            infileS = infileS1;
            outfileS = outfileS1;
            chrRht = chrID4b;
            pos1Rht = pos4b_1;
            pos2Rht = pos4b_2;

        }else if(gene.equals(geneArray[1])){
            infileS = infileS2;
            outfileS = outfileS2;
            chrRht = chrID4d;
            pos1Rht = pos4d_1;
            pos2Rht = pos4d_2;
        }

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            int cnt = 0;
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                if (temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(0));
                int pos = Integer.parseInt(l.get(1));
                if (pos >= (pos1Rht) && pos <= (pos2Rht)){
                    cnt++;
                    bw.write(temp);
                    bw.newLine();
                }
                if (pos > pos2Rht) break;
            }
            System.out.println("Totally " + cnt + " SNPs were disvovered in " + gene);
            br.close();
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void getVCF(){
        String taxaList = "/Users/Aoyue/Documents/TaxaList.txt";


//        String infileS = "/Users/Aoyue/Documents/chr002_S1062.subset.vcf.gz";
//        String outfileS = "/Users/Aoyue/Documents/chr002_spelt.subset.vcf.gz";

        String infileDirS = "/Users/Aoyue/Documents/in";
        String outfileDirS = "/Users/Aoyue/Documents/out";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_subset.vcf").getAbsolutePath();
                CalVCF.extractVCF(infileS,outfileS,taxaList);

                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });


    }


    public void snpAnnotationBuild(){
//        this.mkGeneVCF(); //自己的方法（最终采用）
//        this.mkGeneVCF2(); //来自达兴的方法
//        this.extractInfoFromGeneVCF();
//        this.extractInfoFromGeneVCF_byAoyue();
//        new VariantsSum().mkExonAnnotation2(); //未用
//        new VariantsSum().addAncestral();
//        this.addDAF();
        new VariantsSum().addGerp();
        /**
         * sift 计算
         */
//        new SIFT();
//        this.calFileLine();

//        new VariantsSum().addSift_20210717();
//        new VariantsSum().addDerived_SIFT();
//        new VariantsSum().mergeExonSNPAnnotation();
    }

    /**
     * 获取642倍体中的MAF 个数
     */
    public void getMAFcountfromPop(){

        /**
         * test
         */
//        String infileS = "/Users/Aoyue/Documents/in/chr005.vcf";
//        String outfileS = "/Users/Aoyue/Documents/out/chr005.txt";
//        String taxaInfoDB = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_DD.txt";
//        CalVCF.calMAFcountfromPop(infileS,outfileS,taxaInfoDB);

//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/039_MAFcount/aabbdd";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20210810";
//        String taxaInfoDB = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_AABBDD.txt";
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};


//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/039_MAFcount/aabb";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20210810";
//        String taxaInfoDB = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_AABB.txt";
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/039_MAFcount/dd";
        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20210810";
        String taxaInfoDB = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_DD.txt";
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

//        for (int i = 0; i < chrArr.length; i++) {
//            String chr = chrArr[i];
//            String infileS = new File(infileDirS,"chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chr + "_vmap2.1_MAFcount.txt").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 054_calMAFcountFromPop.jar " +  infileS + " " + outfileS + " " + taxaInfoDB + " > log_20210810_dd.txt");
////            CalVCF.calMAFcountfromPop(infileS,outfileS,taxaInfoDB);
////            System.out.println(chr + "\t finished");
//        }

        SplitScript.splitScript2("/Users/Aoyue/Documents/sh.sh",14,3);
    }




    public void getVCFbyPolyploid(){

        //*** 提取四倍体的VCF ***//
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrAB_vamp2.1_500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/tetra.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABB.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);

        //*** AABBDD - AB sub ***//
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrAB_vmap2.1.500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/hexaploid_ABsub.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABBDD.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);

        //*** DD - D sub ***//
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrD_vmap2.1_500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/Diploid.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_DD.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);

        //*** AABBDD - D sub ***//
//                String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrD_vmap2.1.500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/hexaploid_Dsub.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABBDD.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);


        //*** AABBDD ***//
//        String infileDirS = "/data4/home/aoyue/vmap2/daxing/analysis/007_vmap2_1062/001_subsetVmap2.1/001_vmap2_500k/001_byChrID";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/001_hexaploid";
//        String taxaList = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_AABBDD.txt";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_049_20210808";

        //*** AABBDD ***//
        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/001_outputVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/003_VCF_AABBDD";
        String taxaList = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABBDD.txt";
        String logDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/010_Rht/003_VCF_AABBDD/log";

//
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        String[] chrArr = {"chr021_RhtB1_vmap2.1.vcf.gz","chr023_RhtD1_vmap2.1.vcf.gz"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            String infileS = new File(infileDirS,"chr" + chr + "_vmap2.1.500k.vcf.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chr + "_vmap2.1_hexaploid_subset.vcf.gz").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("nohup java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaList + " > " + logfileS );

            String infileS = new File(infileDirS,chr).getAbsolutePath();
            String outfileS = new File(outfileDirS,chr).getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            CalVCF.extractVCF(infileS,outfileS,taxaList);
        }

//        SplitScript.splitScript2("/Users/Aoyue/Documents/sh.sh",20,3);




    }


    /**
     * 过滤GeneVCF中alt 含有N的位点
     */
    public void filterN_fromVCF(){
//        String infileDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/003_geneSNPVCF_RemoveN";

//        String infileDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0";
//        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/202_VMap2.0";

        String infileDirS = "/data1/home/feilu/filter2";
        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/200_VMap2.0";

//        String infileDirS = "/Users/Aoyue/Documents/in";
//        String outfileDirS = "/Users/Aoyue/Documents/out";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0] + ".vcf").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int cntTotalVariants = 0;
                int cntSNP = 0;
                int cntIndel = 0;
                int cntInsertion = 0;
                int cntDeletion = 0;
                int cntN = 0;
                int cntAllincludeN = 0;

                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) {
                        bw.write(temp);
                        bw.newLine();
                    } else if(temp.startsWith("#C")){
                        List<String> l = PStringUtils.fastSplit(temp);
                        bw.write("#CHROM");
                        StringBuilder sb = new StringBuilder();
                        for (int i = 1; i < l.size(); i++) {
                            sb.append("\t").append(l.get(i));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    else {
                        cntAllincludeN++;
                        List<String> l = PStringUtils.fastSplit(temp);
                        String alt = l.get(4);
                        if (alt.equals("N")) {
                            cntN++;
                            continue; //alt是index为4的列
                        }

                        cntTotalVariants++;

                        if (!alt.equals("D") && !alt.equals("I") && !alt.equals("N")){
                            cntSNP++;
                        }
                        if (alt.equals("D")){
                            cntDeletion++;
                            cntIndel++;
                        }
                        if (alt.equals("I")){
                            cntInsertion++;
                            cntIndel++;
                        }
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName().substring(3,6) + "\t" + cntTotalVariants + "\t" + cntSNP + "\t" + cntIndel + "\t" + cntInsertion + "\t" + cntDeletion + "\t" + cntN);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        // java -jar GeneticLoad.jar > log_removeNfromVCF_20210717.txt 2>&1 &
        // java -jar GeneticLoad.jar > log_removeNfromVMap2.0_20210806.txt 2>&1 &


    }


    public void calFileLine (){
//        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/002_sift/003_output/001_alt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/002_sift/003_output/002_ref";
        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge";

        File[] fsArray = AoFile.getFileArrayInDir(infileDirS);
        for (int i = 0; i < fsArray.length; i++) {
            String filename = fsArray[i].getName();
            String infileS = fsArray[i].getAbsolutePath();
            String chr = filename.substring(3,6);
            int lineNum =  AoFile.countFileRowNumber_withHeader(infileS);
            System.out.println(chr + "\t" + lineNum);
        }
    }

    public void extractVCF(){

        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/008_WheatVMap2_GermplasmInfo_20210708.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subset/chrDsubgenome.100k.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF/diploid.vcf.gz";
        List<String> taxaList = new ArrayList<>();
        RowTable t = new RowTable(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String vcfID = t.getCellAsString(i,0);
            String genomeType = t.getCellAsString(i,4);
            if (genomeType.equals("DD")){
                taxaList.add(vcfID);
            }
        }
        CalVCF.extractVCF(infileS2,outfileS,taxaList);


//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/009_WheatVMap2_GermplasmInfo_20210708.txt";
//        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subset/chrABsubgenome.100k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF/hexaploid.vcf.gz";
//        List<String> taxaList = new ArrayList<>();
//        RowTable t = new RowTable(infileS);
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            String vcfID = t.getCellAsString(i,0);
//            String genomeType = t.getCellAsString(i,4);
//            if (genomeType.equals("AABBDD")){
//                taxaList.add(vcfID);
//            }
//        }
//
//
//        CalVCF.extractVCF(infileS2,outfileS,taxaList);

    }



    /**
     * 老师的方法
     */
    public void addDAF () {
//        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
//        String dirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNPByChr";
        String dirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNP_Annotation";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String header = null;
            List<String> recordList = new ArrayList();
            String tem = null;
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tDAF");
                header = sb.toString();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                }
                br.close();
                BufferedWriter bw = AoFile.writeFile(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                float daf = -1;
                String subMajor = null;
                String subMinor = null;
                String ancestral = null;
                String derivedSIFT = null;
                double subMaf = -1;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t");
                    l = PStringUtils.fastSplit(recordList.get(i));
                    ancestral = l.get(9);
                    if (l.get(5).equals(ancestral)) { //major
                        sb.append((float)Double.parseDouble(l.get(7)));
                    }
                    else if (l.get(6).equals(ancestral)) { //minor
                        sb.append((float)(1- Double.parseDouble(l.get(7))));
                    }
                    else sb.append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                System.out.println(tem);
                e.printStackTrace();
            }
        });
        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addDAFbyFeisMethod_20200720.txt 2>&1 &
        // java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_addDAFbyFeisMethod_20210716.txt 2>&1 &
        // java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_addDAFbyFeisMethod_20210806.txt 2>&1 &

    }

    /**
     * 由于FEI写的程序中，直接从INFO提取 MAF 和判断 Major Minor 会出现错误，故自己写程序计算该列信息
     */
    public void extractInfoFromGeneVCF_byAoyue () {
        /**
         * inputFile
         */
//        String outDirS = "/Users/Aoyue/Documents/out";
//        String vmapDirS = "/Users/Aoyue/Documents/in";
//        String geneHCFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";

        String outDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNP_Annotation";
        String vmapDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz";


        /**
         * parameters need to modify
         */
        File[] fs  = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        AoFile.readheader(geneHCFileS);
        RowTable<String> t = new RowTable<>(geneHCFileS);

        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        int chrNum = chrSet.size(); //chr的个数
        TIntList[] startLists = new TIntList[chrNum]; //list类型的数组，每个数组存放一堆list值
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }

        for (int i = 0; i < t.getRowNumber(); i++) {
            int ifunique = Integer.parseInt(t.getCell(i,3));
            if (ifunique == 0) continue; //过滤不是unique的基因
            startLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 6))); //获取开始位点
            endLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 7))); //获取终止位点
            tranLists[Integer.parseInt(t.getCell(i, 5))-1].add(t.getCell(i, 4)); //获取基因名字
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst("_gene_vmap2.1.vcf", "_geneSNPAnno.txt")).getAbsolutePath();
            StringBuilder sb = new StringBuilder();
            sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tTranscript");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(temp);
                    currentPos = Integer.parseInt(l.get(1));
                    // 开始判断位点是否在基因中
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在开始位点搜索是否在
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;
                    /**
                     * 计算 maf major minor 信息
                     */
                    List<String> lTaxaGeno = new ArrayList<>();
                    String ref = l.get(3);
                    String alt = l.get(4);
                    for (int i = 9; i < l.size(); i++) { //无论有无基因型，都加进去了
                        lTaxaGeno.add(l.get(i));
                    }
                    String[] taxaGenoArray = lTaxaGeno.toArray(new String[lTaxaGeno.size()]);
                    List<String> mafmajorMinorInfo = CalVCF.getMAFMajorMinor(taxaGenoArray,ref,alt);

                    /**
                     * 开始写文件
                     */
                    cnt++;
                    //######## ID ########################### Chr ######################## Pos ####################### Ref
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    //##################### Alt ##################################### Major ######################################### minor ############################## maf 四位小数
                    sb.append("\t").append(l.get(4)).append("\t").append(mafmajorMinorInfo.get(1)).append("\t").append(mafmajorMinorInfo.get(2)).append("\t").append(mafmajorMinorInfo.get(0));
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + "\t" + cnt  + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        //在HPC上运行： nohup java -jar GeneticLoad.jar > log_extractInfoFromVMap2_20210806.txt 2>&1 &
    }


    public void extractInfoFromGeneVCF () {
        /**
         * inputFile
         */
//        String outDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNPByChr";
//        String vmapDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
//        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz";

        String outDirS = "/Users/Aoyue/Documents/out";
        String vmapDirS = "/Users/Aoyue/Documents/in";
        String geneHCFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";

        /**
         * parameters need to modify
         */
        int subLength = 200;
        File[] fs  = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        AoFile.readheader(geneHCFileS);
        RowTable<String> t = new RowTable<>(geneHCFileS);

        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        int chrNum = chrSet.size(); //chr的个数
        TIntList[] startLists = new TIntList[chrNum]; //list类型的数组，每个数组存放一堆list值
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }

        for (int i = 0; i < t.getRowNumber(); i++) {
            int ifunique = Integer.parseInt(t.getCell(i,3));
            if (ifunique == 0) continue; //过滤不是unique的基因
            startLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 6))); //获取开始位点
            endLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 7))); //获取终止位点
            tranLists[Integer.parseInt(t.getCell(i, 5))-1].add(t.getCell(i, 4)); //获取基因名字
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_geneSNP.txt")).getAbsolutePath();
            StringBuilder sb = new StringBuilder();
            sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tTranscript");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) { //如果该行的字符长度小于 200，那么最长设置为该行的实际长度
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    // 开始判断位点是否在基因中
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在开始位点搜索是否在
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;

                    cnt++;
                    //######## ID ########################### Chr ######################## Pos ####################### Ref
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    //##################### Alt ##################
                    sb.append("\t").append(l.get(4)).append("\t");
                    //start - 弃用提取 INFO 的方法
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]);
                    //end - 弃用提取 INFO 的方法

                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + "\t" + cnt  + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        //在HPC上运行： nohup java -jar GeneticLoad.jar > log_extractInfoFromVMap2_20210716.txt 2>&1 &
    }


    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkGeneVCF2 () {
        /**
         * inputFile
         */
        String vmapDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF"; //modify

        RowTableTool<String> rowTable = new RowTableTool<>(hcGeneFileS);
        Predicate<List<String>> p = l->Integer.parseInt(l.get(3))==0;
        rowTable.removeIf(p);

        List<String> geneList = rowTable.getColumn(0);
        PGF pgf = new PGF(geneFeatureFileS);
        Predicate<PGF.Gene> predicate = gene -> !geneList.contains(gene);
        pgf.removeIf(predicate);

        pgf.sortGeneByGeneRange();
        List<File> files = IOTool.getFileListInDirEndsWith(vmapDirS, ".vcf.gz");
        String[] outNames=  files.stream().map(File::getName).map(s -> s.replaceAll("_vmap2.0.vcf.gz","_gene_vmap2.0.vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outputDirS, outNames[e]))) {
                String line, subLine;
                int chrID; int pos;
                List<String> temp;
                while ((line= br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    subLine=line.substring(0, 100);
                    temp=PStringUtils.fastSplit(subLine);
                    chrID =Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    int geneIndex=pgf.getGeneIndex(chrID, pos);
                    if (geneIndex < 0) continue;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_mkGeneVCF_20210715.txt 2>&1 &
    }

    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkGeneVCF () {
        /**
         * inputFile
         */
        String vmapDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF"; //modify
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();

        RowTable<String> t = new RowTable<>(hcGeneFileS);
        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        List<Integer> chrList = new ArrayList<>();
        for (int i = 0; i < chrSet.size(); i++) {
            chrList.add(i+1);
        }
        chrList.parallelStream().forEach(chrID -> {
            String inputVCF = new File(vmapDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_vmap2.0.vcf.gz").getAbsolutePath();
            String outputVCF = new File (outputDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_gene_vmap2.1.vcf").getAbsolutePath();
            List<String> geneList = new ArrayList<>();
            List<String> tranList = new ArrayList<>();
            TIntArrayList tranStartList = new TIntArrayList();
            TIntArrayList tranEndList = new TIntArrayList();
            //获取当前 chr 包含的所有 gene 和 trans
            for (int i = 0; i < t.getRowNumber(); i++) {
                int currentChr = Integer.parseInt(t.getCell(i, 5));
                int ifunique = Integer.parseInt(t.getCell(i,3));
                if (ifunique == 0) continue; //过滤不是unique的基因
                if (currentChr < chrID) continue;
                else if (currentChr > chrID) break;
                geneList.add(t.getCell(i, 0));
                tranList.add(t.getCell(i, 4));
                tranStartList.add(t.getCellAsInteger(i,6));
                tranEndList.add(t.getCellAsInteger(i,7));
            }

            //获取每个基因对应的最长转录本的起始终止位点，然后加入大库的 all gene list
            int[] starts = tranStartList.toArray();
            int[] ends = tranEndList.toArray();

            try {
                BufferedReader br = AoFile.readFile(inputVCF);
                BufferedWriter bw = AoFile.writeFile(outputVCF);
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")) {
                    bw.write(temp); bw.newLine();
                }
                bw.write(temp); bw.newLine(); //#CHROM 这一行
                List<String> l = new ArrayList<>();
                int index = -1;
                int pos = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 100));
                    pos = Integer.parseInt(l.get(1));
                    String alt = l.get(4);
                    if (alt.equals("D") || alt.equals("I") || alt.equals("N")) continue;
                    index = Arrays.binarySearch(starts, pos);
                    if (index < 0) index = -index - 2;
                    if (index < 0) continue;
                    if (pos < ends[index]) {
                        cnt++;
                        bw.write(temp);bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(chrID + "\t" + cnt  + "\t" + "mkGeneVCF");

            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_mkGeneVCF_20210715.txt 2>&1 &
        //java -Xms50g -Xmx200g -jar GeneticLoad_2.jar > log_mkGeneVCF_20210806.txt 2>&1 &
    }

    public void rename(){
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("mv chr" + chr + ".vmap2.vcf.gz chr" + chr + "_vmap2.0.vcf.gz");
            System.out.println("mv chr" + chr + ".vmap2.vcf.gz.tbi chr" + chr + "_vmap2.0.vcf.gz.tbi");
        }
    }


    /**
     * 构建 VMap2.1 使 alt 保持为1，只含有 A T G C 这几种类型
     */
    public void filterSNPtoBi_parallel() {
        String infileDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0";
        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
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
        System.out.println("Chr\tBiallelicNum\tInsertionNum\tDelectionNum\tTri-alleleNum\tN_Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                String outfileS = new File(outfileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int biallelicNum = 0;
                int deletion = 0;
                int insertion = 0;
                int triNum = 0;
                int NANum =0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        String alt = PStringUtils.fastSplit(temp).get(4);
                        if (alt.contains(",")) {
                            triNum++;
                            continue;
                        }
                        if (alt.equals("D")) {
                            deletion++;
                            continue;
                        }
                        if (alt.equals("I")) {
                            insertion++;
                            continue;
                        }
                        if(alt.equals("N")) {
                            NANum++;
                            continue;
                        }

                        biallelicNum++;
                        bw.write(temp);
                        bw.newLine();
                    }
                }

                br.close();
                bw.flush();
                bw.close();
                System.out.println(chrint + "\t" + biallelicNum + "\t" +  insertion + "\t" + deletion + "\t" + triNum + "\t" + NANum);
//                System.out.println(chrint + "\t" + biallelicNum + "\t" + f.getName() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }


    /**
     * 对Fei过滤生成的filter2 进行解压缩，用bgzip压缩，再建立索引
     */
    public void bgzip(){
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("nohup gunzip chr" + chr + ".vmap2.vcf.gz 2>&1 &");
//            System.out.println("nohup bgzip -@ 20 chr" + chr + ".vmap2.vcf && tabix -p vcf chr" + chr + ".vmap2.vcf.gz &");
//            System.out.println("nohup bgzip chr" + chr + "_vmap2.0.vcf && tabix -p vcf chr" + chr + "_vmap2.0.vcf.gz &");
//            System.out.println("nohup tabix -p vcf chr" + chr + "_vmap2.0.vcf.gz &");
            System.out.println("nohup tabix -p vcf chr" + chr + ".vmap2.vcf.gz &");

//            System.out.println("nohup bgzip chr" + chr + ".vmap2.vcf 2>&1 &");
//            System.out.println("nohup bgzip chr" + chr + "_vmap2.1.vcf 2>&1 &");

        }
    }

}
