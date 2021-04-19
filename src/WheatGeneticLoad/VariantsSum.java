/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.*;
import analysis.wheat.VMap2.VMapDBUtils;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.*;
import gnu.trove.set.hash.TIntHashSet;
import nonapi.io.github.classgraph.utils.FileUtils;
import org.apache.commons.math3.stat.inference.TestUtils;
import pgl.graph.tSaw.TablesawUtils;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import pgl.infra.window.SimpleWindow;
import tech.tablesaw.api.IntColumn;
import com.google.common.collect.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class VariantsSum {

    /**
     * 2020年又重新做了一批VMap2,即根据老师写的PopDepth程序计算的每个位点的深度，然后进行70%的密度过滤，取六四二倍体的交集点作为Pos库，
     */

    public void variantsSumFromRebuildVCF(){
//        this.extractInfoFromVMap2();
//        this.mkExonVCF();
//        this.mkExonAnnotation(); //弃用
//        this.mkExonAnnotation2();
//        this.addSift(); //分为2种策略，其中一种已经注释掉。 只添加 alt, 以及同时添加 alt 和 ref 的结果
//        this.addAncestral();
//        this.addDAF_parallel();
//        this.addGerp();
//        this.mergeExonSNPAnnotation();

        //********************** for calculation ****************//
//        this.getDAFtable();
//        this.statisticCodingSNP();
//        this.statisticNonsynSNP();
//        this.getDeleteriouscount();
//        this.getDeleteriousAnnotation();
//        this.countDeleteriousSNP_bySub();
//        this.getGERPdistrbutionFile();
//        this.addGroupToExonAnnotation();

        //*********** count variants in genes *****************//
//        this.countVariantsinGene();
//        this.sortAndFilter();
        //*********** gene summary ********************//
//        this.test();

//*********** 根据 ref sift score 和 alt score 再重新建一遍库 ********************//
// *********** 根据 ref sift score 和 alt score 再重新建一遍库 ********************//
// *********** 根据 ref sift score 和 alt score 再重新建一遍库 ********************//

//        this.mkExonAnnotation2();
//        this.addSift(); //分为2种策略，其中一种已经注释掉。 只添加 alt, 以及同时添加 alt 和 ref 的结果
//        this.addAncestral();
//        this.addDerived_SIFT();
//        this.addDAF_parallel(); //********* 方法有误，不做分析。**********
//        this.addDAF();
//        this.addDAF222222222();  //
//        this.remove01();
//        this.addGerp();
//        this.addRecombination();
//        this.mergeExonSNPAnnotation();

        //********************** for calculation ****************//
//        this.getDAFtable();
//        this.statisticCodingSNP();
//        this.statisticNonsynSNP();
//        this.getDeleteriouscount();
//        this.getDeleteriousAnnotation();
//        this.countDeleteriousSNP_bySub();
//        this.getGERPdistrbutionFile();
//        this.addGroupToExonAnnotation();

        //*********** count variants in genes *****************//
//        this.countVariantsinGene();
//        this.sortAndFilter();

        //*********** ratio of del/syn and nonsyn/syn on genome landscape ********************//
//        this.WindowDelvsSyn_fromExonAnnotation();
//        this.WindowDel_Nonsyn_vsSyn_fromExonAnnotation();
//        this.addRecombinationfromScience(); //该方法凑效！思路：将 del nonysn syn 数据的滑窗设置成和science一致，然后再将重组率文件合并，后续进行其他处理。本次数据分析采用此方法。

//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/003/chr006_SNP_anno.txt");

    }



    /**
     * 老师的方法,copy 一遍
     * 该方法的输入文件是 popdepth 过滤前（2020年5月份之前的结果），chr001有45113sites的exonAnnotation文件，
     * 今日2020年12月30日，重新计算DAF值
     * 为验证popDepth过滤后，Ref bias 有所减轻
     */
    public void addDAF222222222 () {
//        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/001";
//        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/003_popDepthBP_addDAF";

        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002_popDepthBP_addDerivedSIFT_updateAAF";

        List<File> fList = AoFile.getFileListInDir(dirS);
        fList.parallelStream().forEach(f -> {
            String header = null;
            List<String> recordList = new ArrayList();
            String tem = null;
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tDAF\tDAF_ABD\t").append(l.get(9).replaceFirst("AAF", "DAF"));
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
                float dafABD = -1;
                float dafOther = -1;
                String subMajor = null;
                String subMinor = null;
                String ancestral = null;
                String derivedSIFT = null;
                double subMaf = -1;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t");
                    l = PStringUtils.fastSplit(recordList.get(i));
                    ancestral = l.get(15);

                    if (l.get(5).equals(ancestral)) { //major
                        sb.append((float)Double.parseDouble(l.get(7))).append("\t");
                    }
                    else if (l.get(6).equals(ancestral)) { //minor
                        sb.append((float)(1- Double.parseDouble(l.get(7)))).append("\t");
                    }
                    else sb.append("NA\t");


                    if (!l.get(8).equals("NA")){
                        if (Double.parseDouble(l.get(8)) < 0.5) { //AAF_ABD 小于0.5， alt是minor  //******** 总结：要看在亚群体内，谁是major,谁是minor. 与大群体没有关系。
                            subMajor = l.get(3);
                            subMinor = l.get(4);
                            subMaf = Double.parseDouble(l.get(8));
                        }
                        else { //AAF_ABD 大于0.5， alt是major
                            subMajor = l.get(4);
                            subMinor = l.get(3);
                            subMaf = 1 - Double.parseDouble(l.get(8));
                        }
                        if (ancestral.equals(subMajor)) { //alt是minor，ref是major, derived 是 aaf;   alt 是major,ref 是minor, derived是 1-aaf
                            sb.append((float)subMaf).append("\t");
                        }
                        else if (ancestral.equals(subMinor)) { //alt是minor,ref是major,derived是1-aaf; alt是major,ref是mionr,derived是aaf
                            sb.append((float)(1-subMaf)).append("\t");
                        }
                        else sb.append("NA\t");

                    }else{
                        sb.append("NA\t");
                    }


                    tem = recordList.get(i);
                    if (!l.get(9).equals("NA")){
                        if (Double.parseDouble(l.get(9)) < 0.5) {
                            subMajor = l.get(3);
                            subMinor = l.get(4);
                            subMaf = Double.parseDouble(l.get(9));
                        }
                        else {
                            subMajor = l.get(4);
                            subMinor = l.get(3);
                            subMaf = 1 - Double.parseDouble(l.get(9));
                        }
                        if (ancestral.equals(subMajor)) {
                            sb.append((float)subMaf);
                        }
                        else if (ancestral.equals(subMinor)) {
                            sb.append((float)(1-subMaf));
                        }
                        else sb.append("NA");

                    }else{
                        sb.append("NA");
                    }


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
    }



    public void addRecombinationfromScience(){
        String inputFile1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/003_del_nonsyn_syn/001_del_nonsyn_synOnChr_10000000Window1000000step_addEffectiveCDSLength.txt"; //
        String inputFile2 = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_recombination_rate_addScalePos.txt"; //重组率文件
        String outFile = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/003_del_nonsyn_syn/002_del_nonsyn_synOnChr_10000000Window1000000step.txt"; //
        Table<String,String,String> table2=RowTableTool.getTable(inputFile2, 0, 1, 4);
        try (BufferedReader br1 = IOTool.getReader(inputFile1);
             BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            String line=br1.readLine();
            bw.write(line+"\tRecombinationRate");
            bw.newLine();
            List<String> temp;
            String chr, pos;
            while ((line=br1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chr="chr"+temp.get(0);
                pos=temp.get(1);
                if (table2.contains(chr,pos)){
                    temp.add(table2.get(chr,pos));
                }else {
                    temp.add("NA");
                }
                bw.write(String.join("\t", temp));
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 将 del nonsyn syn 的结果写在同一个文件中，不进行分开， 和addCDSLengthInWindow方法的区别在于，前者只计算了del 和 syn 的count 和cds
     * @param dbFileS
     * @param windowSize
     * @param windowStep
     * @param outfile2S
     */
    public void addCDSLengthInWindow2 (String dbFileS, int windowSize, int windowStep, String outfile2S) {
//        int windowSize = 2000000;
//        int windowStep = 1000000;

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String dbFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_2000000Window1000000step.txt";
        String hcGeneFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        RowTable<String> gt = new RowTable<>(hcGeneFileS);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName(); //通过名字排序
        List<String> chromosomeList = RefV1Utils.getChromosomeList(); //返回1A 2A - 7D
        TIntArrayList cdsWindowList = new TIntArrayList(); //

        for (int i = 0; i < chromosomeList.size(); i++) {
            int chrlength = RefV1Utils.getChromosomeLength(chromosomeList.get(i));
            SimpleWindow sw = new SimpleWindow(chrlength, windowSize, windowStep); //new 一个 SimpleWindow 类
            int chrID = RefV1Utils.getChrID(chromosomeList.get(i), 1); //根据chromosomeList中的1A等，获取chrID
            int geneIndex = -1;
            int tranIndex = -1;
            int cdsStart = -1;
            int cdsEnd = -1;
            for (int j = 0; j < gt.getRowNumber(); j++) { //gt已通过名字进行了排序
                int currentChrID = Integer.parseInt(gt.getCell(j, 2));
                if (currentChrID < chrID) continue; //如果，gt文件中的chrID小于循环内的chromosomeList.get(i)，则跳出循环，一直到gt列表里的chrID等于此时正在循环的染色体号
                else if (currentChrID > (chrID +1)) break; //如果，gt文件中的currentChrID 等于chrID+1，说明还是在1A内，其他情况都终止循环
                geneIndex = gf.getGeneIndex(gt.getCell(j, 0)); //根据基因名字（不是转录本的名字，没有.后缀12）获取索引
                for (int k = 0; k < gf.getTranscriptNumber(geneIndex); k++) { //获取最长转录本的index
                    if (!gt.getCell(j,1).equals(gf.getTranscriptName(geneIndex, k)))continue; //如果 gt.getCell(j,1) 是该基因的最长转录本，已经提前总结出
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
            bw.write(header + "\tCDSLength\tDelCountPerSite\tNonsynCountPerSite\tSynCountPerSite");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int line = 0;
            int cdsLength = -1;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(temp);
                cdsLength = cdsWindowList.get(line);
                line++;
                sb.append(temp).append("\t").append(cdsLength).append("\t");
                if (cdsLength != 0) {
                    sb.append((float)((double)Integer.parseInt(l.get(4))/cdsLength)).append("\t");
                    sb.append((float)((double)Integer.parseInt(l.get(5))/cdsLength)).append("\t");
                    sb.append((float)((double)Integer.parseInt(l.get(6))/cdsLength));
                }
                else {
                    sb.append(Float.NaN).append("\t").append(Float.NaN).append("\t").append(Float.NaN);
                }
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

    /**
     * 这里将 del 和 nonsyn 和 syn 的个数都计算出来
     */
    public void WindowDel_Nonsyn_vsSyn_fromExonAnnotation(){

//        int windowSize = 10000000; //10 M
//        int windowStep = 1000000; //1 M

        // 为了进行circos画图
        int windowSize = 20000000; //10 M
        int windowStep = 5000000; //1 M


        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/003_del_nonsyn_syn/001_del_nonsyn_synOnChr_" + windowSize + "Window" + windowStep + "step.txt";
        String outfile2S = outfileS.replaceFirst(".txt","_addEffectiveCDSLength.txt");
        AoFile.readheader(infileS);
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int chrNum = chrArr.length;
        int[][] delePos = new int[chrNum][];
        int[][] synPos = new int[chrNum][];
        int[][] nonsynPos = new int[chrNum][];

        TIntArrayList[] delposList = new TIntArrayList[chrNum];
        TIntArrayList[] synposList = new TIntArrayList[chrNum];
        TIntArrayList[] nonsynposList = new TIntArrayList[chrNum];

        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            delposList[i] = new TIntArrayList();
            synposList[i] = new TIntArrayList();
            nonsynposList[i] = new TIntArrayList();
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
                String variantType = l.get(12);
                String sift = l.get(16); // derived_sift
                String gerp = l.get(20);
                if (variantType.equals("NONSYNONYMOUS")){
                    nonsynposList[index].add(posonchromosome);

                    if (!sift.startsWith("N")){
                        if(!gerp.startsWith("N")){
                            double gerpd = Double.parseDouble(gerp);
                            double siftd = Double.parseDouble(sift);
                            if (siftd <= 0.05 && gerpd >= 1){
                                delposList[index].add(posonchromosome);
                            }
                        }
                    }
                }
                if (variantType.equals("SYNONYMOUS")){
                    synposList[index].add(posonchromosome);
                }
            }
            System.out.println("======== completing the posList DB on all chromosome.");

            BufferedWriter bw = AoFile.writeFile(outfileS);
            String outheader = "CHROM\tBIN_START\tBIN_END\tBIN_START_scale\tDelCount\tNonsynCount\tSynCount\tDelSynRatio\tNonsynSynRatio"; //和vcftools的格式保持一致
            bw.write(outheader);bw.newLine();

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < chrArr.length; i++) {
                String chromosome = chrArr[i];
                int chrLength = RefV1Utils.getChromosomeLength(chromosome);
                SimpleWindow sw = new SimpleWindow(chrLength, windowSize, windowStep);

                sw.addPositionCount(delposList[i].toArray());
                int[] delWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(synposList[i].toArray());
                int[] synWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(nonsynposList[i].toArray());
                int[] nonsynWindowCount = sw.getWindowValuesInt();

                int[] windowStarts = sw.getWindowStarts();
                int[] windowEnds = sw.getWindowEnds();
                for (int j = 0; j < windowStarts.length; j++) {
                    sb.setLength(0);
                    String posscale = WheatUtils.getScaledPos(chromosome,windowStarts[j]);
                    sb.append(chromosome).append("\t").append(windowStarts[j]).append("\t").append(windowEnds[j]).append("\t").append(posscale).append("\t");
                    sb.append(delWindowCount[j]).append("\t").append(nonsynWindowCount[j]).append("\t").append(synWindowCount[j]).append("\t");
                    sb.append((float)((double)delWindowCount[j]/synWindowCount[j])).append("\t").append((float)((double)nonsynWindowCount[j]/synWindowCount[j]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println("======== " + chromosome + " is completed.");
            }
            bw.flush();bw.close();
//*********************************** add cds length *******************************************
            this.addCDSLengthInWindow2(outfileS, windowSize, windowStep,outfile2S);
            System.out.println("======== " + "add effective CDS length.");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * 将重组率文件和del文件合并起来进行计算
     * 该方法错误，完全删除
     *
     * @deprecated
     */
    public void mergeFileandAddScalePos(){
//        String oriFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_2000000Window1000000step_2.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/003_delVSsynOnChr_2000000Window1000000step_addRecombination.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/004_delVSsynOnChr_2000000Window1000000step.txt"; // bin file

        String oriFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_20000000Window5000000step_2.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/003_delVSsynOnChr_50000000Window5000000step_addRecombination.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/004_delVSsynOnChr_20000000Window5000000step.txt";

        //获取recombination的集合
        List<String> recombi = AoFile.getStringList(infileS,3);
        try {
            BufferedReader br = AoFile.readFile(oriFileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tBIN_START_Scale\tRecombinationRate");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int line = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String chromosome = l.get(0);
                int posOnchromosome = Integer.parseInt(l.get(1));
                String startScalePos = WheatUtils.getScaledPos(chromosome,posOnchromosome);
                bw.write(temp + "\t" + startScalePos + "\t" + recombi.get(line));bw.newLine();
                line++;
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 计算每个Bin里面的重组率，并条件到bin文件中
     * 该方法并不可行，计算出的重组率在图上是错的，该步骤完全删除
     *
     * @deprecated
     */
    public void calRecombinationtoBinTable(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //exon annotation file
//        String outfileS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/003_delVSsynOnChr_2000000Window1000000step_addRecombination.txt";
//        int window = 2000000;
//        int step = 1000000;

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //exon annotation file
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/003_delVSsynOnChr_50000000Window5000000step_addRecombination.txt";
        int window = 20000000;
        int step = 5000000;


        AoFile.readheader(infileS);
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int chrNum = chrArr.length;
        TDoubleArrayList[] posList = new TDoubleArrayList[chrNum];
        TDoubleArrayList[] recombiList = new TDoubleArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            posList[i] = new TDoubleArrayList();
            recombiList[i] = new TDoubleArrayList();
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
                String recobi = l.get(21);
                if (recobi.startsWith("N"))continue;
                double recom = Double.parseDouble(recobi);
                String chromosome = RefV1Utils.getChromosome(chr,1); //获取 RefChr
                int posonchromosome = RefV1Utils.getPosOnChromosome(chr,pos);
                int index = Arrays.binarySearch(chrArr,chromosome); //染色体的索引号
                posList[index].add((double)posonchromosome);
                recombiList[index].add(recom);
            }
            System.out.println("======== completing the posList DB on all chromosome.");

            BufferedWriter bw = AoFile.writeFile(outfileS);
            String outheader = "CHROM\tBIN_START\tBIN_END\tRecombinationRate"; //和vcftools的格式保持一致
            bw.write(outheader);bw.newLine();

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < chrArr.length; i++) {
                String chromosome = chrArr[i];
                int chrLength = RefV1Utils.getChromosomeLength(chromosome);

                List<String>[] out = Bin.windowstep_posAve(posList[i],recombiList[i],chrLength,window,step);
                for (int j = 0; j < out[0].size(); j++) {
                    sb.setLength(0);
                    double start = Double.parseDouble(out[0].get(j));
                    double end = start + window;
                    sb.append(chrArr[i]).append("\t").append(String.format("%.0f",start)).append("\t").append(String.format("%.0f",end)).
                            append("\t").append(out[2].get(j));
                    bw.write(sb.toString());bw.newLine();
                }
                System.out.println("======== " + chromosome + " is completed.");
            }
            bw.flush();bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }





//        String[] chrArr = {"1A", "1B"};
//        for (int i = 0; i < chrArr.length; i++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_Wild_emmer_genoTable.txt.gz").getAbsolutePath();
//            List<String> headerList = AoFile.getheader(infileS);
//            //CHROM	BIN_START	BIN_END	N_VARIANTS	HETEROZYGOSITY
//            //CHROM	POS	PI362699	PI487260	PI94648	W11	W13	W14	W15	W16	W18	W19	W20	W5	W7	W8
//            for (int j = 2; j < headerList.size(); j++) {
//                String taxa = headerList.get(j);
//                String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_Wild_emmer" + taxa + "_RH_" + window + "_" + step + ".txt").getAbsolutePath();
//                Bin.ResidualHeterozygosity(infileS,chrArr[i],1,j,window,step,outfileS);
//            }
//            System.out.println();
//        }

    }

    public void addCDSLengthInWindow (String dbFileS, int windowSize, int windowStep, String outfile2S) {
//        int windowSize = 2000000;
//        int windowStep = 1000000;

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String dbFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_2000000Window1000000step.txt";
        String hcGeneFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        RowTable<String> gt = new RowTable<>(hcGeneFileS);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName(); //通过名字排序
        List<String> chromosomeList = RefV1Utils.getChromosomeList();
        TIntArrayList cdsWindowList = new TIntArrayList();

        for (int i = 0; i < chromosomeList.size(); i++) {
            int chrlength = RefV1Utils.getChromosomeLength(chromosomeList.get(i));
            SimpleWindow sw = new SimpleWindow(chrlength, windowSize, windowStep);
            int chrID = RefV1Utils.getChrID(chromosomeList.get(i), 1);
            int geneIndex = -1;
            int tranIndex = -1;
            int cdsStart = -1;
            int cdsEnd = -1;
            for (int j = 0; j < gt.getRowNumber(); j++) {
                int currentChrID = Integer.parseInt(gt.getCell(j, 2));
                if (currentChrID < chrID) continue;
                else if (currentChrID > (chrID +1)) break;
                geneIndex = gf.getGeneIndex(gt.getCell(j, 0)); //根据基因名字获取索引
                for (int k = 0; k < gf.getTranscriptNumber(geneIndex); k++) { //获取最长转录本的index
                    if (!gt.getCell(j,1).equals(gf.getTranscriptName(geneIndex, k)))continue;
                    tranIndex = k;
                    break;
                }
                List<Range> cdsList = gf.getCDSList(geneIndex, tranIndex);
                for (int k = 0; k < cdsList.size(); k++) {
                    cdsStart = cdsList.get(k).start;
                    cdsEnd = cdsList.get(k).end;
                    if (currentChrID == chrID + 1) {
                        cdsStart = RefV1Utils.getPosOnChromosome(currentChrID, cdsStart);
                        cdsEnd = RefV1Utils.getPosOnChromosome(currentChrID, cdsEnd);
                    }
                    sw.addPositionCountFromRange(cdsStart, cdsEnd);
                }
            }
            cdsWindowList.add(sw.getWindowValuesInt());
        }

        try {
            BufferedReader br = AoFile.readFile(dbFileS);
            BufferedWriter bw = AoFile.writeFile(outfile2S);
            String header = br.readLine();
            bw.write(header + "\tCDSLength\tDelCountPerSite\tSynCountPerSite");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int line = 0;
            int cdsLength = -1;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(temp);
                cdsLength = cdsWindowList.get(line);
                line++;
                sb.append(temp).append("\t").append(cdsLength).append("\t");
                if (cdsLength != 0) {
                    sb.append((float)((double)Integer.parseInt(l.get(4))/cdsLength)).append("\t");
                    sb.append((float)((double)Integer.parseInt(l.get(5))/cdsLength));
                }
                else {
                    sb.append(Float.NaN).append("\t").append(Float.NaN).append("\t").append(Float.NaN);
                }
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


        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_" + windowSize + "Window" + windowStep + "step.txt";
        String outfile2S = outfileS.replaceFirst(".txt","_addEffectiveCDSLength.txt");
        AoFile.readheader(infileS);
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int chrNum = chrArr.length;
        int[][] delePos = new int[chrNum][];
        int[][] synPos = new int[chrNum][];

        TIntArrayList[] delposList = new TIntArrayList[chrNum];
        TIntArrayList[] synposList = new TIntArrayList[chrNum];

        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            delposList[i] = new TIntArrayList();
            synposList[i] = new TIntArrayList();
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
                String variantType = l.get(12);
                String sift = l.get(16); // derived_sift
                String gerp = l.get(20);
                if (variantType.equals("NONSYNONYMOUS")){
                    if (!sift.startsWith("N")){
                        if(!gerp.startsWith("N")){
                            double gerpd = Double.parseDouble(gerp);
                            double siftd = Double.parseDouble(sift);
                            if (siftd <= 0.05 && gerpd >= 1){
                                delposList[index].add(posonchromosome);
                            }
                        }
                    }
                }
                if (variantType.equals("SYNONYMOUS")){
                    synposList[index].add(posonchromosome);
                }
            }
            System.out.println("======== completing the posList DB on all chromosome.");

            BufferedWriter bw = AoFile.writeFile(outfileS);
            String outheader = "CHROM\tBIN_START\tBIN_END\tBIN_START_scale\tDelCount\tSynCount\tDelSynRatio"; //和vcftools的格式保持一致
            bw.write(outheader);bw.newLine();

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < chrArr.length; i++) {
                String chromosome = chrArr[i];
                int chrLength = RefV1Utils.getChromosomeLength(chromosome);
                SimpleWindow sw = new SimpleWindow(chrLength, windowSize, windowStep);
                sw.addPositionCount(delposList[i].toArray());
                int[] delWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();

                sw.addPositionCount(synposList[i].toArray());
                int[] synWindowCount = sw.getWindowValuesInt();

                int[] windowStarts = sw.getWindowStarts();
                int[] windowEnds = sw.getWindowEnds();
                for (int j = 0; j < windowStarts.length; j++) {
                    sb.setLength(0);
                    String posscale = WheatUtils.getScaledPos(chromosome,windowStarts[j]);
                    sb.append(chromosome).append("\t").append(windowStarts[j]).append("\t").append(windowEnds[j]).append("\t").append(posscale).append("\t");
                    sb.append(delWindowCount[j]).append("\t").append(synWindowCount[j]).append("\t").append((float)((double)delWindowCount[j]/synWindowCount[j]));
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
     * @deprecated
     * because I need to split the exon annotation file to 2 directory.
     * New method will be worked below.
     */
    public void WindowDelvsSyn () {
        String delInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/del";
        String synInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/syn";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/005_genomeWindowDB/windowDB.txt";
        String header = "CHROM\tBIN_START\tBIN_END\tDelCount\tSynCount\tDelSynRatio"; //和vcftools的格式保持一致
        List<String> chromList = RefV1Utils.getChromosomeList(); //指1A 1B 1D
        int windowSize = 2000000; //2 M
        int windowStep = 1000000; //1 M
        StringBuilder sb = new StringBuilder();
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(header);
            bw.newLine();
            List<String> l = new ArrayList<>();
            //目的，为了获取每条染色体的del位点和syn位点的集合
            for (int i = 0; i < chromList.size(); i++) { //一次循环：21条染色体循环
                int chrLength = RefV1Utils.getChromosomeLength(chromList.get(i)); //指 1A 1B 1D 染色体的长度
                TIntArrayList delList = new TIntArrayList();
                TIntArrayList synList = new TIntArrayList();
                int chrID = RefV1Utils.getChrID(chromList.get(i), 1); //根据1A 1 获取 1-42号染色体的ID
                for (int j = 0; j < 2; j++) { //二次循环：每条染色体的2个拆分进行循环
                    chrID+=j;
                    String delFileS = new File(delInfoDirS, "chr"+ PStringUtils.getNDigitNumber(3, chrID)+"_SNP_anno.txt.gz").getAbsolutePath();
                    String synFileS = new File (synInfoDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_SNP_anno.txt.gz").getAbsolutePath();

                    BufferedReader br = AoFile.readFile(delFileS);
                    String temp = br.readLine(); // read header
                    int pos = -1;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0, 50)); //只为了获取位置
                        pos = Integer.parseInt(l.get(2));
                        if (j == 1) pos = RefV1Utils.getPosOnChromosome(chrID, pos); //如果是2的话，需要重新更换
                        delList.add(pos);
                    }
                    br.close();
                    br = AoFile.readFile(synFileS); //开始处理同义突变
                    temp = br.readLine(); // read header
                    pos = -1;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0, 50));
                        pos = Integer.parseInt(l.get(2));
                        if (j == 1) pos = RefV1Utils.getPosOnChromosome(chrID, pos); //如果是2的话，需要重新更换
                        synList.add(pos);
                    }
                }


                SimpleWindow sw = new SimpleWindow(chrLength, windowSize, windowStep);
                sw.addPositionCount(delList.toArray());
                int[] delWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();
                sw.addPositionCount(synList.toArray());
                int[] synWindowCount = sw.getWindowValuesInt();
                int[] windowStarts = sw.getWindowStarts();
                int[] windowEnds = sw.getWindowEnds();
                for (int j = 0; j < windowStarts.length; j++) {
                    sb.setLength(0);
                    sb.append(chromList.get(i)).append("\t").append(windowStarts[j]).append("\t").append(windowEnds[j]).append("\t");
                    sb.append(delWindowCount[j]).append("\t").append(synWindowCount[j]).append("\t").append((float)((double)delWindowCount[j]/synWindowCount[j]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(chromList.get(i));
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     * ①将统计的额结果进行排序，按照稀有突变（rare variants）中从小到大的顺序进行排序，
     * ②若稀有变异的个数为0或者普通变异的个数为0，则略去
     */
    public void sortAndFilter(){
        this.sort_step1();
        this.remove_step2();
    }



    /**
     * 找到 rare 是0的基因，略去，不进行画图
     */
    public void remove_step2(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/003_sort";
        String infileDirS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/001_byMafVariantsType";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/004_001filterby003";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f ->{
            String infileS = f.getAbsolutePath();
            String name = "001_gene_variantsCount_byMAF_variantsType_" + f.getName().split("_total_")[1];
            String infileS2 = new File(infileDirS2,name).getAbsolutePath();
            String outfileS = new File(outfileDirS,name).getAbsolutePath();
            List<String> geneList = new ArrayList<>();
            RowTable<String> t = new RowTable<>(infileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                String gene = t.getCell(i,0);
                int cntRareVariants = t.getCellAsInteger(i,3);
                if (cntRareVariants==0)continue;
                geneList.add(gene);
            }
            Collections.sort(geneList);

            t=new RowTable<>(infileS2);
            boolean[] ifout = new boolean[t.getRowNumber()];
            for (int i = 0; i < t.getRowNumber(); i++) {
                String gene = t.getCell(i,0);
                int index = Collections.binarySearch(geneList,gene);
                if (index < 0)continue;
                ifout[i] = true;
            }
            t.writeTextTable(outfileS,IOFileFormat.TextGzip,ifout);
        });
    }

    public  void sort_step1(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/002_total";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/003_sort";
        List<File> fs = AoFile.getFileListInDir(infileDirS);
        fs.parallelStream().forEach(f ->{
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            RowTableTool<String> rowTable=new RowTableTool<>(infileS);
            Comparator<List<String>> comparator=Comparator.comparing(l->Integer.parseInt(l.get(3)));
            rowTable.sortBy(comparator);
            rowTable.write(outfileS,IOFileFormat.TextGzip);
        });
    }


    /**
     * 对非常有害的exon数据库进行每个基因的变异数目计数，分成2个分组，第一个分组是同义非同义，第二个分组是DAF值大于5% 和小于等于 5%
     * Gene\tCommonVariants\tRareVariants\tVariantsType\tSub
     * Gene\tCommonVariants\tRareVariants\tNonsynonymous\tSub
     */
    public void countVariantsinGene(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = null;
        String outfileS2 = null;

        AoFile.readheader(infileS);

        //************** 需要修改 *******************//
        double mafThreshold = 0.05; //定义common和rare突变的界限
        String ploidy = "AABBDD";
//        String ploidy = "AABB";
//        String ploidy = "DD";

        //************** 需要修改 *******************//

        int dafColumnIndex= Integer.MIN_VALUE; //六倍体的daf所在的列
        if (ploidy.equals("AABBDD")){
            dafColumnIndex= 8; //六倍体的daf所在的列
            outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/001_byMafVariantsType/001_gene_variantsCount_byMAF_variantsType_AABBDD.txt.gz";
            outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/002_total/001_gene_variantsCount_total_AABBDD.txt.gz";

        }else if(ploidy.equals("AABB")){
            dafColumnIndex= 9; //六倍体的daf所在的列
            outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/001_byMafVariantsType/001_gene_variantsCount_byMAF_variantsType_AABB.txt.gz";
            outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/002_total/001_gene_variantsCount_total_AABB.txt.gz";

        }else if(ploidy.equals("DD")){
            dafColumnIndex= 9; //六倍体的daf所在的列
            outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/001_byMafVariantsType/001_gene_variantsCount_byMAF_variantsType_DD.txt.gz";
            outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/010_countGenevariants/002_total/001_gene_variantsCount_total_DD.txt.gz";
        }

        //******** 数组大小统计 ********//
        String[] genesArray = AoFile.getStringArraybySet(infileS,10);
        Arrays.sort(genesArray);
        System.out.println("Total " + genesArray.length + " genes have SNPs");
        String[] groupMaf = {"Common","Rare"};
        String[] groupVariantType = {"Deleterious","Nonsynonymous","Synonymous"};

        String[] group = {"Common-Deleterious","Common-Nonsynonymous","Common-Synonymous","Rare-Deleterious","Rare-Nonsynonymous","Rare-Synonymous"};
        int[][] count = new int[group.length][genesArray.length]; //最终统计画图
        int[][] countbyMaf = new int[groupMaf.length][genesArray.length]; //统计
        int[][] countbyVariantsType = new int[groupVariantType.length][genesArray.length]; //统计

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            double siftd = Double.NaN;
            double gerpd = Double.NaN;
            double dafd = Double.NaN;
            double maf = Double.NaN;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(1));
                String sub = RefV1Utils.getSubgenomeFromChrID(chr);
                if (!ploidy.contains(sub))continue;//说明必须含有该亚基因组
                String trans = l.get(10);
                int index = Arrays.binarySearch(genesArray,trans);
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################
                String daf = l.get(dafColumnIndex);
                if (daf.equals("NA"))continue; //计数的前提是有DAF值，没有的位点都不考虑,
                dafd = Double.parseDouble(daf);
                //如何判断是 common snp 还是 rare snp?  就看 maf 是大于0.005 还是小于等于0.005
                if(dafd > 0.5) {
                    maf = 1 - dafd;
                }else if (dafd <= 0.5){
                    maf = dafd;
                }

                if (maf > mafThreshold ){ //说明是 common 的

                    if (variantType.equals("SYNONYMOUS")){
                        count[2][index]++;
                        countbyVariantsType[2][index]++;
                        countbyMaf[0][index]++;
                    }
                    if (variantType.equals("NONSYNONYMOUS")){
                        count[1][index]++;
                        countbyVariantsType[1][index]++;
                        countbyMaf[0][index]++;
                    }
                    if (!sift.startsWith("N")){
                        siftd = Double.parseDouble(sift);
                        if (!gerp.startsWith("N")){
                            gerpd = Double.parseDouble(gerp);
                            if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                                count[0][index]++;
                                countbyVariantsType[0][index]++;
                            }
                        }
                    }

                } //maf > mafThreshold

                else if (maf <= mafThreshold ){ //说明是 rare 的
                    if (variantType.equals("SYNONYMOUS")){
                        count[5][index]++;
                        countbyVariantsType[2][index]++;
                        countbyMaf[1][index]++;
                    }
                    if (variantType.equals("NONSYNONYMOUS")){
                        count[4][index]++;
                        countbyVariantsType[1][index]++;
                        countbyMaf[1][index]++;
                    }

                    if (!sift.startsWith("N")){
                        siftd = Double.parseDouble(sift);
                        if (!gerp.startsWith("N")){
                            gerpd = Double.parseDouble(gerp);
                            if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                                count[3][index]++;
                                countbyVariantsType[0][index]++;
                            }
                        }
                    }


                } //maf > mafThreshold
            }
            br.close();

            bw.write("Gene\tSub\tCommonVariants\tRareVariants\tVariantsType");
            bw.newLine();
            //        String[] group = {"Common-Deleterious","Common-Nonsynonymous","Common-Synonymous","Rare-Deleterious","Rare-Nonsynonymous","Rare-Synonymous"};
            for (int i = 0; i < genesArray.length; i++) {
                String gene = genesArray[i];
                String sub = gene.substring(8,9);
                if (!ploidy.contains(sub))continue;//说明必须含有该亚基因组
                bw.write(gene + "\t" + sub + "\t"  + count[2][i] + "\t" + count[5][i] + "\tSynonymous\n");
                bw.write(gene + "\t" + sub + "\t" + count[1][i] + "\t" + count[4][i] + "\tNonsynonymous\n");
                bw.write(gene + "\t" + sub + "\t" + count[0][i] + "\t" + count[3][i] + "\tDeleterious\n");
            }
            bw.flush();bw.close();

            bw = AoFile.writeFile(outfileS2);
            //        int[][] countbyMaf = new int[groupMaf.length][genesArray.length]; //统计
            //        int[][] countbyVariantsType = new int[groupVariantType.length][genesArray.length]; //统计
            bw.write("Gene\tSub\tCommonVariants\tRareVariants\tDeleterious\tNonsynonymous\tSynonymous\n");
            for (int i = 0; i < genesArray.length; i++) {
                String gene = genesArray[i];
                String sub = gene.substring(8,9);
                if (!ploidy.contains(sub))continue;//说明必须含有该亚基因组
                bw.write(gene + "\t" + sub + "\t" + countbyMaf[0][i] + "\t" + countbyMaf[1][i] + "\t" +
                        countbyVariantsType[0][i]  + "\t" + countbyVariantsType[1][i]  + "\t" + countbyVariantsType[2][i]);
                bw.newLine();
            }
            bw.flush();bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 为画DAF的分布，对数据库进行添加分组 "Deleterious","GERP-deleterious","Nonsynonymous-tolerant","SIFT-deleterious","Synonymous"
     * 并提取重要的信息
     * 基于sift < 0.05 和 gerp > 1且是非同义突变
     */
    public void addGroupToExonAnnotation(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/009_addGrouptoExonAnnotation/001_exonSNP_anno_addGroup.txt.gz";

        AoFile.readheader(infileS);
        String[] variantGroup = {"Deleterious","GERP-deleterious","Nonsynonymous-tolerant","SIFT-deleterious","Synonymous"};
        String[] subArray = {"A","B","D"};
        int[][] count = new int[variantGroup.length][subArray.length];
        Arrays.sort(subArray);
        int[] total = new int[variantGroup.length];
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            bw.write("Chr\tPos\tDAF\tGenomeType\tGroup\tSub");
            bw.write("DAF\tGenomeType\tGroup\tSub");

            bw.newLine();
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            double siftd = Double.NaN;
            double gerpd = Double.NaN;
            String genomeType = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(1)); //染色体
                String sub = RefV1Utils.getSubgenomeFromChrID(chr);
                if (sub.equals("D")){
                    genomeType = "DD";
                }
                if (!sub.equals("D")){
                    genomeType="AABB";
                }
                int index = Arrays.binarySearch(subArray,sub);
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String dafABD = l.get(16);
                String dafABorD = l.get(17); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18);


                //分情况：sift有值，并且小于0.5
                if (variantType.equals("SYNONYMOUS")){
                    count[4][index]++;
                    total[4]++;
//                    bw.write(chr + "\t" + pos + "\t" + dafABD + "\t" + "AABBDD" + "\t" + "Synonymous" + "\t" + sub );
//                    bw.newLine();
//                    bw.write(chr + "\t" + pos + "\t" + dafABorD + "\t" + genomeType + "\t" + "Synonymous" + "\t" + sub );
//                    bw.newLine();

                    bw.write(dafABD + "\t" + "AABBDD" + "\t" + "Synonymous" + "\t" + sub );
                    bw.newLine();
                    bw.write(dafABorD + "\t" + genomeType + "\t" + "Synonymous" + "\t" + sub );
                    bw.newLine();

                }

                if (!sift.startsWith("N")){
                    siftd = Double.parseDouble(sift);
                    if (siftd >=0.05 && variantType.equals("NONSYNONYMOUS")) {
                        count[2][index]++;
                        total[2]++;
//                        bw.write(chr + "\t" + pos + "\t" + dafABD + "\t" + "AABBDD" + "\t" + "Nonsynonymous-tolerant" + "\t" + sub );
//                        bw.newLine();
//                        bw.write(chr + "\t" + pos + "\t" + dafABorD + "\t" + genomeType + "\t" + "Nonsynonymous-tolerant" + "\t" + sub );
//                        bw.newLine();
                        bw.write(dafABD + "\t" + "AABBDD" + "\t" + "Nonsynonymous-tolerant" + "\t" + sub );
                        bw.newLine();
                        bw.write(dafABorD + "\t" + genomeType + "\t" + "Nonsynonymous-tolerant" + "\t" + sub );
                        bw.newLine();
                    }
                    if (siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                        count[3][index]++;
                        total[3]++;
//                        bw.write(chr + "\t" + pos + "\t" + dafABD + "\t" + "AABBDD" + "\t" + "SIFT-deleterious" + "\t" + sub );
//                        bw.newLine();
//                        bw.write(chr + "\t" + pos + "\t" + dafABorD + "\t" + genomeType + "\t" + "SIFT-deleterious" + "\t" + sub );
//                        bw.newLine();

                        bw.write(dafABD + "\t" + "AABBDD" + "\t" + "SIFT-deleterious" + "\t" + sub );
                        bw.newLine();
                        bw.write(dafABorD + "\t" + genomeType + "\t" + "SIFT-deleterious" + "\t" + sub );
                        bw.newLine();
                    }
                    if (!gerp.startsWith("N")){
                        gerpd = Double.parseDouble(gerp);
                        if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                            count[0][index]++;
                            total[0]++;
//                            bw.write(chr + "\t" + pos + "\t" + dafABD + "\t" + "AABBDD" + "\t" + "Deleterious" + "\t" + sub );
//                            bw.newLine();
//                            bw.write(chr + "\t" + pos + "\t" + dafABorD + "\t" + genomeType + "\t" + "Deleterious" + "\t" + sub );
//                            bw.newLine();
                            bw.write( dafABD + "\t" + "AABBDD" + "\t" + "Deleterious" + "\t" + sub );
                            bw.newLine();
                            bw.write( dafABorD + "\t" + genomeType + "\t" + "Deleterious" + "\t" + sub );
                            bw.newLine();
                        }
                    }
                }
                if (!gerp.startsWith("N")){
                    gerpd = Double.parseDouble(gerp);
                    if (gerpd > 1 && variantType.equals("NONSYNONYMOUS")) {
                        count[1][index]++;
                        total[1]++;
//                        bw.write(chr + "\t" + pos + "\t" + dafABD + "\t" + "AABBDD" + "\t" + "GERP-deleterious" + "\t" + sub );
//                        bw.newLine();
//                        bw.write(chr + "\t" + pos + "\t" + dafABorD + "\t" + genomeType + "\t" + "GERP-deleterious" + "\t" + sub );
//                        bw.newLine();
                        bw.write(dafABD + "\t" + "AABBDD" + "\t" + "GERP-deleterious" + "\t" + sub );
                        bw.newLine();
                        bw.write( dafABorD + "\t" + genomeType + "\t" + "GERP-deleterious" + "\t" + sub );
                        bw.newLine();

                    }
                }
            }
            br.close();
            bw.flush();bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * 为画出GERP的分布，将文件添加分组并提取特定列数
     * 根据SIFT值，在文件最后再添加一列分组信息 Synonymous
     * Nonsynonymous_tolerent Deleterious
     *
     */
    public void getGERPdistrbutionFile() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/008_gerpDistribution/001_exon_gerpValue_addGroup.txt";
        AoFile.readheader(infileS);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write("###The group is based on sift value");bw.newLine();
            bw.write("Chr\tPos\tGerp\tGroup\tSub");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList();
            double siftd = Double.NaN;
            double gerpd = Double.NaN;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(1)); //染色体
                String sub = RefV1Utils.getSubgenomeFromChrID(chr);
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################


                if (variantType.equals("SYNONYMOUS")){
                    bw.write(chr + "\t" + pos + "\t" + gerp + "\t" + "Synonymous" + "\t" + sub );
                    bw.newLine();
                }
                if (!sift.startsWith("N")){
                    siftd = Double.parseDouble(sift);
                    if (siftd >=0.05 && variantType.equals("NONSYNONYMOUS")) {
                        bw.write(chr + "\t" + pos + "\t" + gerp + "\t" + "Nonsynonymous-tolerant"+ "\t" + sub  );
                        bw.newLine();
                    }
                    if (siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                        bw.write(chr + "\t" + pos + "\t" + gerp + "\t" + "Deleterious"+ "\t" + sub  );
                        bw.newLine();
                    }
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
     * 基于sift < 0.05 和 gerp > 1且是非同义突变，以亚基因组为单位，进行分组计数，后续画柱形图。
     * 再输出一个总体的统计
     */
    public void countDeleteriousSNP_bySub_byMAF(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/007_countVariantType/001_variantGroup_Count_bySub.txt";
        String outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/007_countVariantType/001_variantGroup_Count.txt";

        String[] variantGroup = {"Deleterious","GERP-deleterious","Nonsynonymous-tolerant","SIFT-deleterious","Synonymous"};
        String[] subArray = {"A","B","D"};
        int[][] count = new int[variantGroup.length][subArray.length];
        int[] total = new int[variantGroup.length];

        Arrays.sort(subArray);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            double siftd = Double.NaN;
            double gerpd = Double.NaN;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(1)); //染色体
                String sub = RefV1Utils.getSubgenomeFromChrID(chr);
                int index = Arrays.binarySearch(subArray,sub);
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################

                //分情况：sift有值，并且小于0.5
                if (variantType.equals("SYNONYMOUS")){
                    count[4][index]++;
                    total[4]++;
                }
                if (!sift.startsWith("N")){
                    siftd = Double.parseDouble(sift);
                    if (siftd >=0.05 && variantType.equals("NONSYNONYMOUS")) {
                        count[2][index]++;
                        total[2]++;
                    }
                    if (siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                        count[3][index]++;
                        total[3]++;
                    }
                    if (!gerp.startsWith("N")){
                        gerpd = Double.parseDouble(gerp);
                        if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                            count[0][index]++;
                            total[0]++;
                        }
                    }
                }
                if (!gerp.startsWith("N")){
                    gerpd = Double.parseDouble(gerp);
                    if (gerpd > 1 && variantType.equals("NONSYNONYMOUS")) {
                        count[1][index]++;
                        total[1]++;
                    }
                }
            }
            br.close();

            bw.write("VariantGroup\tSub\tCount");
            bw.newLine();
            for (int i = 0; i < variantGroup.length; i++) {
                String group1 = variantGroup[i];
                for (int j = 0; j < subArray.length; j++) {
                    String sub = subArray[j];
                    int num = count[i][j];
                    bw.write(group1 + "\t" + sub + "\t" + num);
                    bw.newLine();
                    System.out.println(group1 + "\t" + sub + "\t" + num);
                }
            }
            bw.flush();bw.close();

            bw = AoFile.writeFile(outfileS2);
            bw.write("VariantGroup\tCount");bw.newLine();
            for (int i = 0; i < variantGroup.length; i++) {
                String group = variantGroup[i];
                int num = total[i];
                bw.write(group + "\t" + num);
                bw.newLine();
                System.out.println(group + "\t" + num);
            }
            bw.flush();bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 基于sift < 0.05 和 gerp > 1且是非同义突变，以亚基因组为单位，进行分组计数，后续画柱形图。
     * 再输出一个总体的统计
     */
    public void countDeleteriousSNP_bySub(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/007_countVariantType/001_variantGroup_Count_bySub.txt";
        String outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/007_countVariantType/001_variantGroup_Count.txt";

        String[] variantGroup = {"Deleterious","GERP-deleterious","Nonsynonymous-tolerant","SIFT-deleterious","Synonymous"};
        String[] subArray = {"A","B","D"};
        int[][] count = new int[variantGroup.length][subArray.length];
        int[] total = new int[variantGroup.length];

        Arrays.sort(subArray);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            double siftd = Double.NaN;
            double gerpd = Double.NaN;

            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(1)); //染色体
                String sub = RefV1Utils.getSubgenomeFromChrID(chr);
                int index = Arrays.binarySearch(subArray,sub);
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################

                //分情况：sift有值，并且小于0.5
                if (variantType.equals("SYNONYMOUS")){
                    count[4][index]++;
                    total[4]++;
                }
                if (!sift.startsWith("N")){
                    siftd = Double.parseDouble(sift);
                    if (siftd >=0.05 && variantType.equals("NONSYNONYMOUS")) {
                        count[2][index]++;
                        total[2]++;
                    }
                    if (siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                        count[3][index]++;
                        total[3]++;
                    }
                    if (!gerp.startsWith("N")){
                        gerpd = Double.parseDouble(gerp);
                        if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                            count[0][index]++;
                            total[0]++;
                        }
                    }
                }
                if (!gerp.startsWith("N")){
                    gerpd = Double.parseDouble(gerp);
                    if (gerpd > 1 && variantType.equals("NONSYNONYMOUS")) {
                        count[1][index]++;
                        total[1]++;
                    }
                }
            }
            br.close();

            bw.write("VariantGroup\tSub\tCount");
            bw.newLine();
            for (int i = 0; i < variantGroup.length; i++) {
                String group1 = variantGroup[i];
                for (int j = 0; j < subArray.length; j++) {
                    String sub = subArray[j];
                    int num = count[i][j];
                    bw.write(group1 + "\t" + sub + "\t" + num);
                    bw.newLine();
                    System.out.println(group1 + "\t" + sub + "\t" + num);
                }
            }
            bw.flush();bw.close();

            bw = AoFile.writeFile(outfileS2);
            bw.write("VariantGroup\tCount");bw.newLine();
            for (int i = 0; i < variantGroup.length; i++) {
                String group = variantGroup[i];
                int num = total[i];
                bw.write(group + "\t" + num);
                bw.newLine();
                System.out.println(group + "\t" + num);
            }
            bw.flush();bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     *  获取满足sift < 0.05, gerp >1 且是 nonsynonymous 的位点的注释信息,并添加一列标注亚基因组的信息
     */
    public void getDeleteriousAnnotation(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/006_exonAnnnotation_sift0.05Gerp1/001_exonSNP_sift0.05_GERP3_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/006_exonAnnnotation_sift0.05Gerp1/001_exonSNP_sift0.05_GERP1_anno.txt.gz";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header+ "\tSub"); bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(1)); //染色体
                String sub = RefV1Utils.getSubgenomeFromChrID(chr);
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################

                //分情况：sift有值，并且小于0.5
                double siftd = Double.NaN;
                double gerpd = Double.NaN;
                if (!sift.startsWith("N")){
                    siftd = Double.parseDouble(sift);
                    if (!gerp.startsWith("N")){
                        gerpd = Double.parseDouble(gerp);
//                        if (gerpd > 3 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
//                            bw.write(temp + "\t" + sub);
//                            bw.newLine();
//                            cnt++;
//                        }
                        if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) {
                            bw.write(temp + "\t" + sub);
                            bw.newLine();
                            cnt++;
                        }
                    }
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(cnt + "\tdel SNPs ");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 计数
     * ①non-syn 且 sift < 0.05;
     * ②non-syn 且 gerp >1;
     * ③non-syn && sift < 0.05 && gerp >1
     */
    public void getDeleteriouscount(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        String[] variantGroup = {"Nonsyn-sift","Nonsyn-gerp","Nonsyn-SiftGerp"};
        int[] count = new int[variantGroup.length];

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int index = Integer.parseInt(l.get(1)) - 1; //染色体号的索引 ################ 需要修改 需要修改 需要修改 ################
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################

                //分情况：sift有值，并且小于0.5
                double siftd = Double.NaN;
                double gerpd = Double.NaN;
                if (!sift.startsWith("N")){
                    siftd = Double.parseDouble(sift);
                    if (siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) count[0]++;
                    if (!gerp.startsWith("N")){
                        gerpd = Double.parseDouble(gerp);
                        if (gerpd > 1 && siftd < 0.05 && variantType.equals("NONSYNONYMOUS")) count[2]++;
                    }
                }
                if (!gerp.startsWith("N")){
                    gerpd = Double.parseDouble(gerp);
                    if (gerpd > 1 && variantType.equals("NONSYNONYMOUS")) count[1]++;

                }
            }
            br.close();
            for (int i = 0; i < variantGroup.length; i++) {
                System.out.println(variantGroup[i] + "\t" + count[i]);
            }
            //Nonsyn-sift	129135
            //Nonsyn-gerp	131456
            //Nonsyn-SiftGerp	61151

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 统计编码区的SNP个数，因外显子包含5'UTR CDS 3'UTR，我们只统计CDS区域的SNP个数
     */
    public void statisticNonsynSNP(){
        //从总的合并的注释文件中获取信息
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        AoFile.readheader(infileS);
        AoMath.countCaseInGroup(infileS,12);
        //NA	3
        //STOP-LOSS	1138
        //START-LOST	710
        //SYNONYMOUS	381280
        //NONCODING	408589
        //NONSYNONYMOUS	421555
        //STOP-GAIN	9194
    }

    /**
     * 统计编码区的SNP个数，因外显子包含5'UTR CDS 3'UTR，我们只统计CDS区域的SNP个数
     */
    public void statisticCodingSNP(){
        //从总的合并的注释文件中获取信息
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        AoFile.readheader(infileS);
        AoMath.countCaseInGroup(infileS,11);
        //CDS	813880
        //UTR_3	329116
        //UTR_5	79473
    }


    /**
     * 根据计算的DAF值，进行分bin,得到表格进行画图操作
     * 根据数据库动态创建分组，将该分组内的所有数字建立list，进行bin的统计，并返回每个bin的比例
     *
     */
    public void getDAFtable() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014/002_basedGerpPhyloP";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014/003_basedSIFT_ratio";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/018_getDAFtablefrom014/005_basedonlyGERP";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/012_exonSNPAnnotation_merge_filterHeter0.05";

//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/004_DAFtable"; //总共的ABD
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/005_DATtable_barley_urartu";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/004_DAFtable_barley_secale_parsimony";

//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/006_DAFtable_barley_secale_parsimony_filterHeter0.05";

        //********************  需要手动设置 START ****************************//
//        String infileDirS = ""; 输入文件是合并的EXON Annotation数据库
//        String outfileDirS = ""; //输出文件的目录，即会有6个文件 A sub B sub  D sub  only ABD only AB only D
//        new File(outfileDirS).mkdirs();
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz");

        //********************  需要手动设置 START ****************************//

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/"; //001_exonSNP_anno.txt.gz
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/005_DAFtable_barley_secale_parsimony";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/005_DAFtable_barley_secale_parsimony/001_DAFtable_siftGERP";
//        new File(outfileDirS).mkdirs();
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz");

        ////*************** 用 derivedSIFT 重建数据库后的重计算
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge"; //输入文件是合并的EXON Annotation数据库
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/005_DAFtable_barley_secale_parsimony/002_DAFtable_basedSiftGerp_newDBuseDerivedSIFT"; //输出文件的目录，即会有6个文件 A sub B sub  D sub  only ABD only AB only D
        new File(outfileDirS).mkdirs();
        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation/chr001_SNP_anno.txt.gz");


        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            try {
                //************************************ 第一阶段，定义输出输出文件，读写文件 ************************//
                String infileS = f.getAbsolutePath();
                //********************  需要手动设置 START ****************************//
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_onlyABD.txt").getAbsolutePath(); //只能画总体的ABD六倍体
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_onlyAB.txt").getAbsolutePath(); //只能画总体的AB四倍体
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_onlyD.txt").getAbsolutePath(); //只能画总体的D二倍体
//
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_Asubgenome.txt").getAbsolutePath(); //只有A亚基因组的结果
//                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_Bsubgenome.txt").getAbsolutePath(); //只有B亚基因组的结果
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "binTable_Dsubgenome.txt").getAbsolutePath(); //只有D亚基因组的结果
                //********************  需要手动设置 END ****************************//

                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                //************************************ 第二阶段，创建相关变量，并读入文件判断 ************************//
                String[] group = {"Deleterious SNPs", "Nonsynonymous-tolerant SNPs", "Synonymous SNPs"};
                TDoubleArrayList[] dafABD = new TDoubleArrayList[group.length];
                for (int i = 0; i < dafABD.length; i++) { //初始化List
                    dafABD[i] = new TDoubleArrayList();
                }

                TDoubleArrayList[] dafAB = new TDoubleArrayList[group.length];
                for (int i = 0; i < dafAB.length; i++) { //初始化List
                    dafAB[i] = new TDoubleArrayList();
                }

                TDoubleArrayList[] daf = new TDoubleArrayList[group.length];
                for (int i = 0; i < daf.length; i++) { //初始化List
                    daf[i] = new TDoubleArrayList();
                }

                double sift = Double.NaN;
                double gerp = Double.NaN;
                String temp = null;
                String header = br.readLine();
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) { //我想得到 DAF_ABD的pos集合， DAF_AB的pos集合，以及合并数据的DAF集合。每个集合又分为3类，一类是同义突变，一类是非同义突变，sift值小于0.05，一类是非同义突变，sfft大于0.05
                    l = PStringUtils.fastSplit(temp);
                    int chrID = Integer.parseInt(l.get(1));
                    int posID = Integer.parseInt(l.get(2));
                    String chr = RefV1Utils.getChromosome(chrID,posID); //根据 chr pos 获取该位点坐在的亚基因组

                    //********************  需要手动设置 START ****************************//
//                    if(chr.contains("D"))continue; //只能用于AABB总体画图
//                    if(chr.contains("A") || chr.contains("B"))continue; //只能用于DD总体画图

//                    if(chr.contains("B") || chr.contains("D"))continue; //只能用于A亚基因组
//                    if(chr.contains("A") || chr.contains("D"))continue; //只能用于B亚基因组
                    if(chr.contains("A") || chr.contains("B"))continue; //只能用于D亚基因组
                    //********************  需要手动设置 END ****************************//

                    //********************  需要手动设置 START ****************************//
                    System.out.println(temp);
                    String type = l.get(12);
//                    String siftscore = l.get(13);
//                    String gerpscore = l.get(18);
//                    String DAF_ABD = l.get(16); //大麦黑麦Parsimony法为ancestral allele计算的DAF值
//                    String DAF_AB = l.get(17); //大麦黑麦Parsimony法为ancestral allele计算的DAF值
//                    String DAF = l.get(15); //大麦黑麦Parsimony法为ancestral allele计算的DAF值

                    String siftscore = l.get(16);
                    String gerpscore = l.get(20);
                    String DAF_ABD = l.get(18); //大麦黑麦Parsimony法为ancestral allele计算的DAF值
                    String DAF_AB = l.get(19); //大麦黑麦Parsimony法为ancestral allele计算的DAF值
                    String DAF = l.get(17); //大麦黑麦Parsimony法为ancestral allele计算的DAF值
//********************  需要手动设置 END ****************************//
                    //如果变异类型是同义突变，那么就不用做任何判断；直接加上分组 Synonymous 并写入
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值小于0.5，gerp和phylop存在，且gerp大于1，且phylop大于0.5；那么加上分组 Deleterious 并写入； gerp 值和 phylop值不满足条件的，那么就不进行分组
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值大于0.5，那么加上分组 Nonsynonymous_tolerent 并写入
                    //如果变异类型是非同义突变，SIFT值不存在，那么不分组 不写入
                    if (type.equals("SYNONYMOUS")) {
                        if (!DAF.startsWith("N")) { //说明是有值的
                            double value = Double.parseDouble(DAF);
                            daf[2].add(value);
                        }
                        if (!DAF_ABD.startsWith("N")) { //说明是有值的
                            double value = Double.parseDouble(DAF_ABD);
                            dafABD[2].add(value);
                        }
                        if (!DAF_AB.startsWith("N")) { //说明是有值的
                            double value = Double.parseDouble(DAF_AB);
                            dafAB[2].add(value);
                        }

                    }
                    if (type.equals("NONSYNONYMOUS")) {

//                        if (!gerpscore.startsWith("N")) { //均有值存在
//                            gerp = Double.parseDouble(gerpscore);
//                            if (gerp > 1) {
//                                if (!DAF.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF);
//                                    daf[0].add(value);
//                                }
//                                if (!DAF_ABD.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_ABD);
//                                    dafABD[0].add(value);
//                                }
//                                if (!DAF_AB.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_AB);
//                                    dafAB[0].add(value);
//                                }
//                            }else if (gerp <= 1){ //说明是 非同义突变并且gerp小于1
//                                if (!DAF.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF);
//                                    daf[1].add(value);
//                                }
//                                if (!DAF_ABD.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_ABD);
//                                    dafABD[1].add(value);
//                                }
//                                if (!DAF_AB.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_AB);
//                                    dafAB[1].add(value);
//                                }
//
//                            }
//                        }

                        if (!siftscore.startsWith("N")) {
                            sift = Double.parseDouble(siftscore);
                            if (sift <= 0.05) {
                                //添加gerp phyloP分组信息 情况一：
//                                if (!gerpscore.startsWith("N") && (!phylopscore.startsWith("N"))) { //均有值存在
//                                    gerp = Double.parseDouble(gerpscore);
//                                    phylop = Double.parseDouble(phylopscore);
//                                    if (gerp > 1 && (phylop > 0.5)) {
//                                        if (!DAF.startsWith("N")) { //说明是有值的
//                                            double value = Double.parseDouble(DAF);
//                                            daf[0].add(value);
//                                        }
//                                        if (!DAF_ABD.startsWith("N")) { //说明是有值的
//                                            double value = Double.parseDouble(DAF_ABD);
//                                            dafABD[0].add(value);
//                                        }
//                                        if (!DAF_AB.startsWith("N")) { //说明是有值的
//                                            double value = Double.parseDouble(DAF_AB);
//                                            dafAB[0].add(value);
//                                        }
//                                    }
//                                }


                                //不添加gerp phyloP分组信息 情况二：
//**************************** 可供选择 *********************************************** //
//                                if (!DAF.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF);
//                                    daf[0].add(value);
//                                }
//                                if (!DAF_ABD.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_ABD);
//                                    dafABD[0].add(value);
//                                }
//                                if (!DAF_AB.startsWith("N")) { //说明是有值的
//                                    double value = Double.parseDouble(DAF_AB);
//                                    dafAB[0].add(value);

//**************************** 可供选择 *********************************************** //

                                //添加gerp 分组信息 情况三：
                                if (!gerpscore.startsWith("N")) { //均有值存在
                                    gerp = Double.parseDouble(gerpscore);
                                    if (gerp > 1) {
                                        if (!DAF.startsWith("N")) { //说明是有值的
                                            double value = Double.parseDouble(DAF);
                                            daf[0].add(value);
                                        }
                                        if (!DAF_ABD.startsWith("N")) { //说明是有值的
                                            double value = Double.parseDouble(DAF_ABD);
                                            dafABD[0].add(value);
                                        }
                                        if (!DAF_AB.startsWith("N")) { //说明是有值的
                                            double value = Double.parseDouble(DAF_AB);
                                            dafAB[0].add(value);
                                        }
                                    }
                                }

                            } else { //sift值大于0.05

                                if (!DAF.startsWith("N")) { //说明是有值的
                                    double value = Double.parseDouble(DAF);
                                    daf[1].add(value);
                                }
                                if (!DAF_ABD.startsWith("N")) { //说明是有值的
                                    double value = Double.parseDouble(DAF_ABD);
                                    dafABD[1].add(value);
                                }
                                if (!DAF_AB.startsWith("N")) { //说明是有值的
                                    double value = Double.parseDouble(DAF_AB);
                                    dafAB[1].add(value);
                                }
                            } //if (sift <= 0.05) {
                        } //if (!siftscore.startsWith("N")) {
                    } //if (type.equals("NONSYNONYMOUS")) {
                }
                br.close(); //文件阅读完毕！！

                //************************************ 第三阶段，开始返回bin的比例 ************************//

                List[] l1 = this.mkBarplotofDAF(dafABD[0], 20); //六倍体的有害突变
                List[] l2 = this.mkBarplotofDAF(dafABD[1], 20); //六倍体的可忍受非同义突变
                List[] l3 = this.mkBarplotofDAF(dafABD[2], 20); //六倍体的同义突变
                List[] l4 = this.mkBarplotofDAF(dafAB[0], 20); //四倍体的有害突变
                List[] l5 = this.mkBarplotofDAF(dafAB[1], 20); //四倍体的可忍受非同义突变
                List[] l6 = this.mkBarplotofDAF(dafAB[2], 20); //四倍体的同义突变
                List[] l7 = this.mkBarplotofDAF(daf[0], 20); //二倍体的有害突变
                List[] l8 = this.mkBarplotofDAF(daf[1], 20); //二倍体的可忍受非同义突变
                List[] l9 = this.mkBarplotofDAF(daf[2], 20); //二倍体的同义突变
                //************************************ 第四阶段，开始写出文件 ************************//
                bw.write("Xaxes\tDAF_ABD\tDAF_AB\tDAF\tGroup");
                bw.newLine();
                for (int i = 0; i < l1[0].size(); i++) { //本行都是有害突变
                    bw.write(String.valueOf(l1[0].get(i)) + "\t" + String.valueOf(l1[1].get(i)) + "\t" + String.valueOf(l4[1].get(i))+ "\t" + String.valueOf(l7[1].get(i)) + "\t" + group[0]);
                    bw.newLine();
                }
                for (int i = 0; i < l1[0].size(); i++) { //本行都是可忍受非同义突变
                    bw.write(String.valueOf(l1[0].get(i))+ "\t" + String.valueOf(l2[1].get(i)) + "\t" + String.valueOf(l5[1].get(i))+ "\t" + String.valueOf(l8[1].get(i))+ "\t" + group[1]);
                    bw.newLine();
                }
                for (int i = 0; i < l1[0].size(); i++) { //本行都是同义突变
                    bw.write(String.valueOf(l1[0].get(i)) + "\t"+ String.valueOf(l3[1].get(i)) + "\t" + String.valueOf(l6[1].get(i))+ "\t" + String.valueOf(l9[1].get(i))+ "\t" + group[2]);
                    bw.newLine();
                }


//                bw.write("Xaxes\tDeleterious_SNPs\tNonsynonymous_tolerant_SNPs\tSynonymous_SNPs\tGroup");
//                bw.newLine();
//                for (int i = 0; i < l1[0].size(); i++) {
//                    bw.write(String.valueOf(l1[0].get(i)) + "\t" + String.valueOf(l1[1].get(i)) + "\t" + String.valueOf(l2[1].get(i))+ "\t" + String.valueOf(l3[1].get(i)) + "\t" + "Hexaploid");
//                    bw.newLine();
//                }
//                for (int i = 0; i < l1[0].size(); i++) {
//                    bw.write(String.valueOf(l1[0].get(i))+ "\t" + String.valueOf(l4[1].get(i)) + "\t" + String.valueOf(l5[1].get(i))+ "\t" + String.valueOf(l6[1].get(i))+ "\t" + "Tetraploid");
//                    bw.newLine();
//                }
//                for (int i = 0; i < l1[0].size(); i++) {
//                    bw.write(String.valueOf(l1[0].get(i)) + "\t"+ String.valueOf(l7[1].get(i)) + "\t" + String.valueOf(l8[1].get(i))+ "\t" + String.valueOf(l9[1].get(i))+ "\t" + "Diploid");
//                    bw.newLine();
//                }

                bw.flush();
                bw.close();

                //
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    public List[] mkBarplotofDAF(TDoubleArrayList dafList, int bins) {
        List[] l = new List[2];
        for (int i = 0; i < l.length; i++) {
            l[i] = new ArrayList();
        }
//        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble("1");

        //先建立bound数组,假如 bin=20,那么每个bin的大小是0.05 bound[0] =0; bound[1]=0.05, bound[19]=0.95 即只取了区间的左边开始位置
        double[] bound = new double[bins];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double) length / bins * i;
        }

        double[] daf = new double[bins]; //确定落入每个区间的个数，进行计算，20个bin，20个count数字
        for (int i = 0; i < dafList.size(); i++) { //对每个daf值进行区间的判断
            double value = dafList.get(i);
            int index = Arrays.binarySearch(bound, value);
            if (index < 0) {
                index = -index - 2;
            }
            //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
            //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
            daf[index]++; //值落入第i种变异的第index个区间的个数
        }
        //开始计算每个区间落入点的比例
        for (int i = 0; i < daf.length; i++) {
            daf[i] = daf[i] / dafList.size();
        }
        System.out.println("This list is " + dafList.size() + "  size");
        //开始写出文件
        try {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < bound.length; i++) {
                String coordinate = String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2);
                String proportion = String.format("%.4f", daf[i]);
                l[0].add(coordinate);
                l[1].add(proportion);
                System.out.println(coordinate + "   " + proportion);
            }


        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return l;
    }

    public void mergeExonSNPAnnotation(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/014_exonSNPAnnotation";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        AoFile.mergeTxt(infileDirS,outfileS);

//        String infileDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
//        String outfileS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        AoFile.mergeTxt(infileDirS,outfileS);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        AoFile.mergeTxt(infileDirS,outfileS);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/003";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002_exonSNP_anno_beforePopDepth.txt.gz";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/004_popDepthBP_remove01";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002_exonSNP_anno_beforePopDepth.txt.gz";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/004_popDepthBP_addDerivedSIFT_updateAAF_addDAF_remove01";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002_exonSNP_anno_beforePopDepth.txt.gz";

        AoFile.mergeTxt(infileDirS,outfileS);


        //java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_mergeExonSNPAnnotation_20200609.txt 2>&1 &
    }

    /**
     * 为单个位点添加重组率的结果
     */
    public void addRecombination () {
        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/012/001_exonSNPAnnotation";
        String outDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/012/002";
        String recombinationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_recombination_rate_chrID.txt";
        ColumnTable<String> t = new ColumnTable<>(recombinationFileS);
        int chrNum = Integer.parseInt(t.getCell(t.getRowNumber()-1, 0)); //合计42条染色体
        TIntArrayList[] startLists = new TIntArrayList[chrNum];
        TIntArrayList[] endLists = new TIntArrayList[chrNum];
        TFloatArrayList[] crossLists = new TFloatArrayList[chrNum];
        for (int i = 0; i < startLists.length; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            crossLists[i] = new TFloatArrayList();
        }
        int index = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            index = Integer.parseInt(t.getCell(i, 0))-1;
            startLists[index].add(Integer.parseInt(t.getCell(i, 1)));
            endLists[index].add(Integer.parseInt(t.getCell(i, 2)));
            crossLists[index].add(Float.parseFloat(t.getCell(i, 3)));
        }
        List<File> fList = IOUtils.getFileListInDirEndsWith(dirS, ".gz");
        fList.parallelStream().forEach(f -> {
            //String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> recordList = two.getSecondElement();
            String[] tem = header.split("\t");
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                StringBuilder sb = new StringBuilder(header);
                sb.append("\t").append("RecombinationRate");
                bw.write(sb.toString());
                bw.newLine();
                int chrIndex = -1;
                int posIndex = -1;
                int currentPos = -1;
                List<String> l  = null;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(recordList.get(i));
                    chrIndex = Integer.parseInt(l.get(1))-1;
                    currentPos = Integer.parseInt(l.get(2));
                    posIndex = startLists[chrIndex].binarySearch(currentPos);
                    if (posIndex < 0) posIndex = -posIndex-2;
                    if (posIndex < 0) {
                        sb.append(recordList.get(i)).append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    if (currentPos < endLists[chrIndex].get(posIndex)) {
                        sb.append(recordList.get(i)).append("\t").append(crossLists[chrIndex].get(posIndex));
                    }
                    else {
                        sb.append(recordList.get(i)).append("\t").append("NA");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void addGerp () {
        String gerpDirS = "/data4/home/aoyue/vmap2/feilu/003_annotation/003_gerp/byChr_29way";

//        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";

//        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";

        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/004_exonAnnotation_beforeDepthFilter/002";
        List<File> fList = AoFile.getFileListInDir(dirS);
        fList.parallelStream().forEach(f -> {
            String gerpFileS = f.getName().split("_")[0]+"_gerp.txt.gz";
            gerpFileS = new File (gerpDirS, gerpFileS).getAbsolutePath();
            String header = null;
            List<String> recordList = new ArrayList();
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                header = br.readLine();
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tGerp");
                header = sb.toString();
                String temp = null;
                TIntArrayList posList = new TIntArrayList();
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                    l = PStringUtils.fastSplit(temp);
                    posList.add(Integer.parseInt(l.get(2)));
                }
                br.close();
                br = AoFile.readFile(gerpFileS);
                br.readLine(); //header
                int pos = -1;
                int index = -1;
                String[] gerp = new String[posList.size()];
                for (int i = 0; i < gerp.length; i++) gerp[i] = "NA";
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    pos = Integer.parseInt(l.get(1));
                    index = posList.binarySearch(pos);
                    if (index < 0) continue;
                    gerp[index] = l.get(2);
                }
                br.close();
                BufferedWriter bw = AoFile.writeFile(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                for (int i = 0; i < posList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t").append(gerp[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(f.getAbsolutePath() + " is completed.");
        });
        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addGerp_20200608.txt 2>&1 &
        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addGerp_20200721.txt 2>&1 &

    }

    /**
     * 在上文fei计算的DAF基础上，将DAF中等于1或者0的值都替换成NA，不做分析。
     */
    public void remove01(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/012/001_exonAnnotation";
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/012/002_remove01";
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/012/001_exonAnnotation/chr001_SNP_anno.txt.gz");

        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/003_popDepthBP_addDerivedSIFT_updateAAF_addDAF/chr001_SNP_anno_updateAAF.txt");
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/003_popDepthBP_addDAF";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/004_popDepthBP_remove01";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/003_popDepthBP_addDerivedSIFT_updateAAF_addDAF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/004_popDepthBP_addDerivedSIFT_updateAAF_addDAF_remove01";


        List<File> fsList = AoFile.getFileListInDir(infileDirS);

        fsList.parallelStream().forEach(f ->{
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            try{
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                bw.write(header);bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while((temp=br.readLine()) != null){
                    l = PStringUtils.fastSplit(temp);
//                    String DAF = l.get(17);
//                    String DAF_ABD = l.get(18);
//                    String DAF_AB = l.get(19);

                    String DAF = l.get(18);
                    String DAF_ABD = l.get(19);
                    String DAF_AB = l.get(20);
                    String daf = AoMath.replace10toNA(DAF);
                    String dafABD = AoMath.replace10toNA(DAF_ABD);
                    String dafAB = AoMath.replace10toNA(DAF_AB);
                    //####### 是NA的保持不变，是 1 或者 0 的统统替换成 NA 写个方法把~
                    for (int i = 0; i < 18; i++) {
                        bw.write(l.get(i) + "\t");
                    }
                    bw.write(daf + "\t" + dafABD + "\t" + dafAB);
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 老师的方法
     */
    public void addDAF () {
        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
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
                sb.append("\tDAF\tDAF_ABD\t").append(l.get(9).replaceFirst("AAF", "DAF"));
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
                float dafABD = -1;
                float dafOther = -1;
                String subMajor = null;
                String subMinor = null;
                String ancestral = null;
                String derivedSIFT = null;
                double subMaf = -1;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t");
                    l = PStringUtils.fastSplit(recordList.get(i));
                    ancestral = l.get(15);

                    if (l.get(5).equals(ancestral)) { //major
                        sb.append((float)Double.parseDouble(l.get(7))).append("\t");
                    }
                    else if (l.get(6).equals(ancestral)) { //minor
                        sb.append((float)(1- Double.parseDouble(l.get(7)))).append("\t");
                    }
                    else sb.append("NA\t");
                    if (Double.parseDouble(l.get(8)) < 0.5) { //AAF_ABD 小于0.5， alt是minor  //******** 总结：要看在亚群体内，谁是major,谁是minor. 与大群体没有关系。
                        subMajor = l.get(3);
                        subMinor = l.get(4);
                        subMaf = Double.parseDouble(l.get(8));
                    }
                    else { //AAF_ABD 小于0.5， alt是major
                        subMajor = l.get(4);
                        subMinor = l.get(3);
                        subMaf = 1 - Double.parseDouble(l.get(8));
                    }
                    if (ancestral.equals(subMajor)) { //alt是minor，ref是major, derived 是 aaf;   alt 是major,ref 是minor, derived是 1-aaf
                        sb.append((float)subMaf).append("\t");
                    }
                    else if (ancestral.equals(subMinor)) { //alt是minor,ref是major,derived是1-aaf; alt是major,ref是mionr,derived是aaf
                        sb.append((float)(1-subMaf)).append("\t");
                    }
                    else sb.append("NA\t");

                    tem = recordList.get(i);
                    if (Double.parseDouble(l.get(9)) < 0.5) {
                        subMajor = l.get(3);
                        subMinor = l.get(4);
                        subMaf = Double.parseDouble(l.get(9));
                    }
                    else {
                        subMajor = l.get(4);
                        subMinor = l.get(3);
                        subMaf = 1 - Double.parseDouble(l.get(9));
                    }
                    if (ancestral.equals(subMajor)) {
                        sb.append((float)subMaf);
                    }
                    else if (ancestral.equals(subMinor)) {
                        sb.append((float)(1-subMaf));
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
    }


    public void addDAF_parallel(){ //本地运行常用
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_exonSNPAnnotation";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/004_exonSNPAnnotation";

        String infileDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/004";

        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            this.calDAF(f.getAbsolutePath(),outfileS);
            System.out.println(f.getName() + "\tis completed at " + outfileS);
        });

        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addFDA_20200720.txt 2>&1 &
    }

    /**
     * Goal:根据 ancestral allele，计算Daf,Daf_ABD Daf_AB Daf_D
     * ************ 注意注意注意：该方法存在理解性错误！！！！！！ 该方法错误错误错误错误！！！！！ 请勿使用！！！！ *********
     */
    public void calDAF(String dbfileS, String outfileS) { //String dbfileS, String ancS, String outfileS
//        String dbfileS = "";
//        String outfileS = "";
        boolean ifd = false;
        double daf = Double.NaN;
        double daf_ABD = Double.NaN;
        double daf_AB = Double.NaN;
        int cntAncNum = 0;
        File f = new File(dbfileS); //根据ancestral allele 文件，得到染色体号
        int chr = Integer.parseInt(f.getName().substring(3, 6));
        //根据染色体号进行AB还是D的判断
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, chr) > -1) { //说明是属于D的
            ifd = true;
        }
        try {
            String chrS = PStringUtils.getNDigitNumber(3, chr);
            BufferedReader br = AoFile.readFile(dbfileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine(); //read header
            if (ifd == false) {
                bw.write(temp + "\tDaf_barley_secaleParsimony\tDaf_ABD_barley_secaleParsimony\tDaf_AB_barley_secaleParsimony");
                bw.newLine();
            } else if (ifd == true) {
                bw.write(temp + "\tDaf_barley_secaleParsimony\tDaf_ABD_barley_secaleParsimony\tDaf_D_barley_secaleParsimony");
                bw.newLine();
            }

            int cntAnc = 0;
            int cntAncNotMajororMinor = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(2));
                String major = l.get(5);
                String minor = l.get(6);
                double maf = Double.parseDouble(l.get(7));
                double AAF_ABD = Double.parseDouble(l.get(8));
                double AAF_AB = Double.parseDouble(l.get(9));

                //################### 需要修改 //###################//###################//###################//###################
//               String ancAllele = l.get(22); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancAllele = l.get(15); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
//                String ancAllele = l.get(14);
                //################### 需要修改 //###################//###################//###################//###################

                StringBuilder sb = new StringBuilder();
                if (!ancAllele.equals("NA")) { //表明含有anc
                    //如果ancestral allele存在,且等于major，则derived allele等于minor, daf 就等于maf
                    //如果ancestral allele存在,且等于minor，则derived allele等于major, daf 就等于 1-daf1
                    if (ancAllele.equals(minor)) {
                        cntAnc++;
                        daf = 1 - maf;
                        if (AAF_ABD > 0.5) { //说明AAF_ABD是major， 祖先状态是minor的，所有DAF是major
                            daf_ABD = AAF_ABD;
                        } else if (AAF_ABD < 0.5) { //说明AAF_ABD是minor， 祖先状态是minor的，所有DAF是major
                            daf_ABD = 1 - AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = 1 - AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (ancAllele.equals(major)) {
                        cntAnc++;
                        daf = maf;
                        if (AAF_ABD > 0.5) { //说明AAF_ABD是major， 祖先状态是major的，所有DAF是minor
                            daf_ABD = 1 - AAF_ABD;
                        } else if (AAF_ABD < 0.5) {
                            daf_ABD = AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = 1 - AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (!ancAllele.equals(minor) && (!ancAllele.equals(major))) {
                        sb.append(temp).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                        cntAncNotMajororMinor++;
                    }

                } else { //表明不含anc
                    sb.append(temp).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            double ratio = (double) cntAncNotMajororMinor / (cntAncNotMajororMinor + cntAnc);
            bw.flush();
            bw.close();
            br.close();
            System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
            System.out.println(new File(dbfileS).getName() + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor. The ratio is " + String.format("%.4f", ratio));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 根据alt sift score 和 ref sift score 和 ancestral allele 的状态， 判断 derived allele 的得分， 可能是 alt 也可能是 ref
     */
    public void addDerived_SIFT(){
//        String infileDirS = "/Users/Aoyue/Downloads/003_exonSNPAnnotation";
//        String outfileDirS = "/Users/Aoyue/Downloads/003";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/001_popDepthBP";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002_popDepthBP_addDerivedSIFT";

        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/001_popDepthBP/chr006_SNP_anno.txt");
        System.out.println();
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.forEach(f ->{
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            String header = null;
            try {
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                header = br.readLine();
                bw.write(header + "\tDerived_SIFT");bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String ref = l.get(3);
                    String alt = l.get(4);
                    String ancestral = l.get(15);
                    String refSIFT = l.get(14);
                    String altSIFT = l.get(13);
                    if (ancestral.equals("NA")){
                        bw.write(temp + "\tNA");
                        bw.newLine();
                    }else if (!ancestral.equals("NA")){
                        if (ancestral.equals(ref)){ //说明 derived allele 是 alt
                            bw.write(temp + "\t" + altSIFT);
                            bw.newLine();
                        }
                        if(ancestral.equals(alt)){ //说明 derived allele 是 ref
                            bw.write(temp + "\t" + refSIFT);
                            bw.newLine();
                        }
                        if(!ancestral.equals(ref) && (!ancestral.equals(alt))){ //ancestral 即不等于ref 又不等于 alt
                            bw.write(temp + "\tNA");
                            bw.newLine();
                        }
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addDerivedSIFT_20200720.txt 2>&1 &
    }

    public void addAncestral () {
//        String inDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/002_byChrID";
//        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_exonSNPAnnotation";

        String inDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/003_annotation/001_ancestral";
        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";

        List<File> fList = AoFile.getFileListInDir(inDirS);
        fList.parallelStream().forEach(f -> {
            String annoFileS = f.getName().split("_")[0]+"_SNP_anno.txt";
            annoFileS = new File(dirS, annoFileS).getAbsolutePath();
            String header = null;
            List<String> recordList = new ArrayList();
            TIntArrayList posList = new TIntArrayList();
            String[] ancestral = null;
            try {
                BufferedReader br = AoFile.readFile(annoFileS);
                header = br.readLine();
                String temp = null;
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                    l = PStringUtils.fastSplit(temp);
                    posList.add(Integer.parseInt(l.get(2)));
                }
                ancestral = new String[posList.size()]; //注释库的Pos的信息库
                for (int i = 0; i < ancestral.length; i++) ancestral[i] = "NA";
                br.close();


                br = AoFile.readFile(f.getAbsolutePath()); //read ancestral file
                temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    int index = posList.binarySearch(pos);
                    if (index<0) continue;
                    ancestral[index] = l.get(2);
                }
                br.close();


                BufferedWriter bw = AoFile.writeFile(annoFileS);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tAncestral");
                bw.write(sb.toString());
                bw.newLine();
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t").append(ancestral[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addAncestral_20200720.txt 2>&1 &   cat /data4/home/aoyue/vmap2/aaPlantGenetics/log_addAncestral_20200720.txt
    }

    /**
     * 此方法是同时添加 ref 和 sift 的方法，注意 在ref 的结果中只提取sift score 这一项内容；
     * refsift 和  altsift 的结果，行数和pos的位置都一一对应。
     * 另外注意：单独添加sift alt的方法在下文可以找到，需要用到的时候，请打开注释号。
     * CHROM	POS	REF_ALLELE	ALT_ALLELE	TRANSCRIPT_ID	GENE_ID	GENE_NAME	REGION	VARIANT_TYPE	REF_AMINO	ALT_AMINO	AMINO_POS	SIFT_SCORE	SIFT_MEDIAN	NUM_SEQS	dbSNSIFT_PREDICTION
     * 1	1145382	C	A	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	NONSYNONYMOUS	F	L	68	0.68	2.70	227	novel	TOLERATED
     * 1	1145386	G	T	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	STOP-GAIN	E	*	70	NA	NA	NA	novel	NA
     * 1	1145397	G	A	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	SYNONYMOUS	T	T	73	0.13	2.70	227	novel	TOLERATED
     * 1	1145405	C	T	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	NONSYNONYMOUS	A	V	76	0.03	2.70	227	novel	DELETERIOUS
     * 1	1145406	C	T	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	SYNONYMOUS	A	A	76	1.00	2.70	227	novel	TOLERATED
     * 1	1145418	C	T	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	SYNONYMOUS	S	S	80	1.00	2.70	227	novel	TOLERATED
     * 1	1145442	G	A	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	SYNONYMOUS	K	K	88	1.00	2.71	228	novel	TOLERATED
     * 1	1145448	A	G	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	SYNONYMOUS	Q	Q	90	1.00	2.73	228	novel	TOLERATED
     * 1	1145454	G	C	TraesCS1A02G001800.1	TraesCS1A02G001800	TraesCS1A02G001800	CDS	SYNONYMOUS	V	V	92	1.00	2.73	227	novel	TOLERATED
     */
    public void addSift() {
        String siftAltDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/003_result_Vmap2.1-2020_exonVCF/output"; //chr001_exon_vmap2.1_SIFTannotations.xls.gz
        String siftRefDirS = "/data4/home/aoyue/vmap2/analysis/028_sift/001_exon_refref/output";
        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation"; //库文件
        File[] fs = new File(siftAltDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".xls.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String dbFileS = f.getName().split("_")[0]+"_SNP_anno.txt";
            dbFileS = new File (dirS, dbFileS).getAbsolutePath();
            String refFileS = new File (siftRefDirS, f.getName().split("_")[0]+"_exon_vmap2.1_RefRef_SIFTannotations.xls.gz").getAbsolutePath();
            RowTable<String> tAlt = new RowTable (f.getAbsolutePath());
            RowTable<String> tRef = new RowTable (refFileS);
            SIFTRecord[] records = new SIFTRecord[tAlt.getRowNumber()];
            for (int i = 0; i < records.length; i++) {
                SIFTRecord s = new SIFTRecord(Integer.parseInt(tAlt.getCell(i, 1)), tAlt.getCell(i, 3), tAlt.getCell(i, 4), tAlt.getCell(i,
                        7), tAlt.getCell(i, 8), tAlt.getCell(i, 12), tRef.getCell(i, 12));
                records[i] = s;
            }
            Arrays.sort(records);
            try {
                List<String> dbList = new ArrayList();
                String temp = null;
                BufferedReader br = AoFile.readFile(dbFileS);
                String header = br.readLine();
                while ((temp = br.readLine()) != null) {
                    dbList.add(temp);
                }
                br.close();
                BufferedWriter bw = AoFile.writeFile(dbFileS);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tRegion\tVariant_type\tAlt_SIFT\tRef_SIFT");
                bw.write(sb.toString());
                bw.newLine();
                List<String> l = null;
                for (int i = 0; i < dbList.size(); i++) {
                    l = PStringUtils.fastSplit(dbList.get(i));
                    SIFTRecord query = new SIFTRecord(Integer.parseInt(l.get(2)), l.get(4), l.get(10));
                    int index = Arrays.binarySearch(records, query);
                    if (index < 0) continue;
                    sb.setLength(0);
                    sb.append(dbList.get(i)).append("\t");
                    sb.append(records[index].region).append("\t");
                    sb.append(records[index].type).append("\t");
                    sb.append(records[index].altSift).append("\t").append(records[index].refSift);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addSift_20200720.txt 2>&1 &   cat /data4/home/aoyue/vmap2/aaPlantGenetics/log_addSift_20200720.txt
    }

    class SIFTRecord implements Comparable<SIFTRecord> {

        public int pos;
        public String alt;
        public String transcript;
        public String region;
        public String type;
        public String altSift;
        public String refSift;


        public SIFTRecord(int pos, String alt, String transcript) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
        }

        public SIFTRecord(int pos, String alt, String transcript, String region, String type, String altSift, String refSift) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
            this.region = region;
            this.type = type;
            this.altSift = altSift;
            this.refSift = refSift;
        }

        @Override
        public int compareTo(SIFTRecord o) { //
            if (this.pos < o.pos) {
                return -1;
            } else if (this.pos == o.pos) {
                int index = this.alt.compareTo(o.alt);
                if (index < 0) {
                    return -1;
                }
                else if (index > 0) {
                    return 1;
                }
                else return transcript.compareTo(o.transcript);
            } else {
                return 1;
            }
        }
    }




//    public void addSift() {
//        String siftAltDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/003_result_Vmap2.1-2020_exonVCF/output";
////        String siftRefDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/001_sift/output_ref";
//        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
//        List<File> fList = AoFile.getFileListInDir(siftAltDirS);
//        fList.parallelStream().forEach(f -> {
//            String dbFileS = f.getName().split("_")[0]+"_SNP_anno.txt.gz";
//            dbFileS = new File (dirS, dbFileS).getAbsolutePath();
////            String refFileS = new File (siftRefDirS, f.getName().split("_")[0]+"_exon_vmap2.1_reverseRefAlt_SIFTannotations.xls.gz").getAbsolutePath();
//            RowTable<String> tAlt = new RowTable (f.getAbsolutePath());
////            RowTable<String> tRef = new RowTable (refFileS);
//            SIFTRecord[] records = new SIFTRecord[tAlt.getRowNumber()];
//            for (int i = 0; i < records.length; i++) {
//                SIFTRecord s = new SIFTRecord(Integer.parseInt(tAlt.getCell(i, 1)), tAlt.getCell(i, 3), tAlt.getCell(i, 4), tAlt.getCell(i,
//                        7), tAlt.getCell(i, 8), tAlt.getCell(i, 12));
//                records[i] = s;
//            }
//            Arrays.sort(records);
//            try {
//                List<String> dbList = new ArrayList(); //是原来文件的每一行的结果
//                String temp = null;
//                BufferedReader br = AoFile.readFile(dbFileS);
//                String header = br.readLine(); //
//                while ((temp = br.readLine()) != null) {
//                    dbList.add(temp);
//                }
//                br.close();
//                BufferedWriter bw = AoFile.writeFile(dbFileS);
//                StringBuilder sb = new StringBuilder(header);
//                sb.append("\tRegion\tVariant_type\tAlt_SIFT");
//                bw.write(sb.toString());
//                bw.newLine();
//                List<String> l = null;
//                for (int i = 0; i < dbList.size(); i++) {
//                    l = PStringUtils.fastSplit(dbList.get(i));
//                    SIFTRecord query = new SIFTRecord(Integer.parseInt(l.get(2)), l.get(4), l.get(10));
//                    int index = Arrays.binarySearch(records, query);
//                    if (index < 0) continue;
//                    sb.setLength(0);
//                    sb.append(dbList.get(i)).append("\t");
//                    sb.append(records[index].region).append("\t");
//                    sb.append(records[index].type).append("\t");
//                    sb.append(records[index].altSift);
//                    bw.write(sb.toString());
//                    bw.newLine();
//                }
//                bw.flush();
//                bw.close();
//            }
//            catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addSift_20200608.txt 2>&1 &
//
//
//    }

//    class SIFTRecord implements Comparable<SIFTRecord> {
//
//        public int pos;
//        public String alt;
//        public String transcript;
//        public String region;
//        public String type;
//        public String altSift;
//
//
//
//        public SIFTRecord(int pos, String alt, String transcript) {
//            this.pos = pos;
//            this.alt = alt;
//            this.transcript = transcript;
//        }
//
//        public SIFTRecord(int pos, String alt, String transcript, String region, String type, String altSift) {
//            this.pos = pos;
//            this.alt = alt;
//            this.transcript = transcript;
//            this.region = region;
//            this.type = type;
//            this.altSift = altSift;
//
//        }
//
//        @Override
//        public int compareTo(SIFTRecord o) { //
//            if (this.pos < o.pos) {
//                return -1;
//            } else if (this.pos == o.pos) {
//                int index = this.alt.compareTo(o.alt);
//                if (index < 0) {
//                    return -1;
//                }
//                else if (index > 0) {
//                    return 1;
//                }
//                else return transcript.compareTo(o.transcript);
//            } else {
//                return 1;
//            }
//        }
//    }




    /**
     * 如何根据Pos直接找出对应的基因？
     * 第一步：将pgf文件的基因按照染色体的位置进行排序； 基因排序
     * 第二步：根据chr pos，使用getGeneIndex方法得到基因的index  找基因索引
     * 第三步：使用getGeneName方法，得到基因的名字； 找基因名字
     * 第四步：使用getLongestTranscriptIndex方法，得到该基因最长转录本的index;  找转录本索引
     * 第五步：使用getTranscriptName方法，根据基因的index和最长转录本的index得到转录本的名字；  找转录本名字
     */
    public void mkExonAnnotation2(){
        int subLength = 200; //取VCF文件的前200个字符串
        String vmapDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/002_exonSNPVCF";
        String outDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
        File[] fs = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        //将高置信度的基因文件读进表格
        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/001_geneHC/geneHC.txt";
        RowTable<String> t = new RowTable<>(geneHCFileS);
        //将pgf文件新建对象，并通过基因名字排列
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();

        vmapList.parallelStream().forEach(f -> {
            int chrID = Integer.parseInt(f.getName().substring(3, 6));
            List<String> geneList = new ArrayList<>();
            List<String> tranList = new ArrayList<>();
            TIntArrayList startLists = new TIntArrayList();
            TIntArrayList endLists = new TIntArrayList();


            for (int i = 0; i < t.getRowNumber() ; i++) {
                int currentChr = Integer.parseInt(t.getCell(i, 2));
                if (currentChr < chrID) continue; //因为文件是按照chrpos排列的，所以表格中chr小于当前文件的chr,就继续循环
                else if (currentChr > chrID) break; //表格中chr大于当前的chr，就终止循环，说明已经建立好基因列表
                geneList.add(t.getCell(i, 0));
                tranList.add(t.getCell(i, 1));
                startLists.add(Integer.parseInt(t.getCell(i,3))); //本条染色体内的所有基因的起始位置组合
                endLists.add(Integer.parseInt(t.getCell(i,4)));
            }

            int[] starts = startLists.toArray(new int[startLists.size()]);
            int[] ends = endLists.toArray(new int[endLists.size()]);

            String outfileS = new File(outDirS, f.getName().replaceFirst("_exon_vmap2.1.vcf.gz", "_SNP_anno.txt")).getAbsolutePath();
            int[] dc = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(dc);
            StringBuilder sb = new StringBuilder();
            if (Arrays.binarySearch(dc, chrID) < 0) {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB\tTranscript");
            } else {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D\tTranscript");
            }
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
                int cntkept = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) {
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    posIndex = Arrays.binarySearch(starts,currentPos);
                    if (posIndex < 0) {
                        posIndex = -posIndex - 2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos > ends[posIndex]) continue;

                    cntkept++;
                    String trans = tranList.get(posIndex); //在第几个基因中搜到，就是第几个转录本。 注意：只是本条染色体内部的基因集合
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    sb.append("\t").append(l.get(4)).append("\t");
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""), ",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    } else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]).append("\t").append(ll.get(7).split("=")[1]).append("\t").append(ll.get(8).split("=")[1]);
                    sb.append("\t").append(trans);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed at " + outfileS + " with\t" + cntkept + " SNPs");
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_mkExonAnnotation2_20200607.txt 2>&1 &
        //java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_mkExonAnnotation2_20200720.txt 2>&1 &
        //cat /data4/home/aoyue/vmap2/aaPlantGenetics/log_mkExonAnnotation2_20200720.txt

    }

    /**
     * @deprecated
     */
    public void mkExonAnnotation() {
        String inDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/001_genicSNPByChr";
        String outDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
        String exonVCFDirS = "/data4/home/aoyue/vmap2/feilu/002_genicSNP/002_exonSNPVCF";
        List<File> fList = AoFile.getFileListInDir(inDirS);
        fList.parallelStream().forEach(f -> { //chr001_vmap2.1_genicSNP.txt.gz
            String name = f.getName().substring(3,6);
            int chr = Integer.parseInt(name);
            String outfileS = f.getName().split("_")[0]+"_SNP_anno.txt.gz";
            outfileS = new File (outDirS, outfileS).getAbsolutePath();
            String exonVCFfileS = new File(exonVCFDirS,"chr" + name + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            TIntArrayList posList = CalVCF.extractVCFPos(exonVCFfileS);
            posList.sort();
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                bw.write(header);bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(2));
                    int index = posList.binarySearch(pos);
                    if (index < 0)continue;
                    cnt++;
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + " is completed at " + outfileS + " with\t" + cnt + "SNPs");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            //java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_mkExonAnnotation_20200606.txt 2>&1 &
            //cat /data4/home/aoyue/vmap2/aaPlantGenetics/log_mkExonAnnotation_20200606.txt

        });

    }

    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkExonVCF () {
        String vmapDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/001_geneHC/geneHC.txt"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/002_exonSNPVCF"; //modify
        GeneFeature gf = new GeneFeature(geneFeatureFileS);

        gf.sortGeneByName();
        RowTable<String> t = new RowTable<>(hcGeneFileS);
        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(2)); //get chr的set集合
        List<Integer> chrList = new ArrayList<>();
        for (int i = 0; i < chrSet.size(); i++) {
            chrList.add(i+1);
        }
        chrList.parallelStream().forEach(chrID -> {
            String inputVCF = new File (vmapDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_vmap2.1.vcf").getAbsolutePath();
            String outputVCF = new File (outputDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_exon_vmap2.1.vcf").getAbsolutePath();
            List<String> geneList = new ArrayList<>();
            List<String> tranList = new ArrayList<>();
            //获取当前 chr 包含的所有 gene 和 trans
            for (int i = 0; i < t.getRowNumber(); i++) {
                int currentChr = Integer.parseInt(t.getCell(i, 2));
                if (currentChr < chrID) continue;
                else if (currentChr > chrID) break;
                geneList.add(t.getCell(i, 0));
                tranList.add(t.getCell(i, 1));
            }

            //获取该染色体的每个基因对应的最长转录本的exonlist,然后加入大库的 all exon list。
            int geneIndex = -1;
            List<Range> allexonList = new ArrayList<>();
            for (int i = 0; i < geneList.size(); i++) {
                geneIndex = gf.getGeneIndex(geneList.get(i));
                for (int j = 0; j < gf.getTranscriptNumber(geneIndex); j++) { //该基因的所有转录本的循环
                    if (!tranList.get(i).equals(gf.getTranscriptName(geneIndex, j))) continue; //验证对应的trans是否和pgf文件中的最长转录本一致，不一致，程序退出
                    List<Range> exonList = gf.getExonList(geneIndex, j);
                    allexonList.addAll(exonList);
                }
            }
            Collections.sort(allexonList);
            int[] starts = new int[allexonList.size()];
            int[] ends = new int[allexonList.size()];
            for (int i = 0; i < starts.length; i++) {
                starts[i] = allexonList.get(i).getRangeStart();
                ends[i] = allexonList.get(i).getRangeEnd();
            }
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
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 100));
                    pos = Integer.parseInt(l.get(1));
                    index = Arrays.binarySearch(starts, pos);
                    if (index < 0) index = -index - 2; //zai
                    if (index < 0) continue;
                    if (pos < ends[index]) {
                        bw.write(temp);bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(chrID+"  mkExonVCF");
        });
        //java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_mkExonVCF_20200606.txt 2>&1 &
    }


    /**
     * 提取VMap2的基因区间的基本信息
     */
    public void extractInfoFromVMap2 () {
        int subLength = 200;
        String outDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/001_genicSNPByChr";
        String vmapDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1";
        File[] fs  = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/001_geneHC/geneHC.txt";
        AoFile.readheader(geneHCFileS);
        tech.tablesaw.api.Table t = TablesawUtils.readTsv(geneHCFileS);
        System.out.println(t.structure());
        t.sortAscendingOn("Chr", "TranStart");
        IntColumn chrColumn = t.intColumn("chr");
        int chrNum = chrColumn.countUnique(); //chr的个数
        TIntList[] startLists = new TIntList[chrNum]; //list类型的数组，每个数组存放一堆list值
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }

        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3)));
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4)));
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1));
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_genicSNP.txt")).getAbsolutePath();
            int[] dc = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(dc);
            StringBuilder sb = new StringBuilder();
            if (Arrays.binarySearch(dc, chrIndex+1) < 0) {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB\tTranscript");
            }
            else {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D\tTranscript");
            }
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
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) {
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    posIndex = startLists[chrIndex].binarySearch(currentPos);
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    sb.append("\t").append(l.get(4)).append("\t");
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]).append("\t").append(ll.get(7).split("=")[1]).append("\t").append(ll.get(8).split("=")[1]);
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        //在HPC上运行： java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_extractInfoFromVMap2_20200606.txt 2>&1 &
    }

    public VariantsSum() {

        this.variantsSumFromRebuildVCF();

//        new CountSites().mergeChr1Aand2A_bysubgenome();

//        this.mkSNPsummary("/data4/home/aoyue/vmap2/genotype/mergedVCF/002_biMAF0.005VCF/", "/data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/");
        //this.addAncestralAllele("/Users/Aoyue/Documents/out", "/Users/Aoyue/Documents/out1", "/Users/Aoyue/Documents/out2");
        //this.scriptAddAncAllele();
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/delSNP");
//new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/nonsyTolerantSNP");
//new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/synSNP");
        // new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/004_merge/delSNP");
        // new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/004_merge/nonsyTolerantSNP");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_merge006/004_merge/synSNP");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/delSNP/chrA.subgenome.delSNP.changeChrPos.txt");
        //new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/delSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/delSNP/chrB.subgenome.delSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/nonsyTolerantSNP/chrA.subgenome.nonsyTolerantSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/nonsyTolerantSNP/chrB.subgenome.nonsyTolerantSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/synSNP/chrA.subgenome.synSNP.changeChrPos.txt");
//new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/003_merge1A/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/synSNP/chrB.subgenome.synSNP.changeChrPos.txt");
//this.mkBarplotOfSNPs();

        //this.classifySNPs("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/005_addAncestralAllele/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/");
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/delSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/delSNP");
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/nonsyTolerantSNP");
        //new CountSites().changechrPos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/006_ori/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/synSNP");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/delSNP/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/delSNP/");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/nonsyTolerantSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/nonsyTolerantSNP");
        //new CountSites().mergefile1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/007_changeChrPos/synSNP", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/synSNP");
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/delSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/delSNP/chrD.subgenome.delSNP.changeChrPos.txt");
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/nonsyTolerantSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/nonsyTolerantSNP/chrD.subgenome.nonsyTolerantSNP.changeChrPos.txt");
//    new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/008_merge1A/synSNP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/synSNP/chrD.subgenome.synSNP.changeChrPos.txt");
//        this.mkBarplotOfSNPs();
//        this.mkSNPsummary_step1("/Users/Aoyue/Documents/chr001_vmap2.1_line50.vcf.gz", "/Users/Aoyue/Documents/out/chr001_vmap2.1_line50.annotation.txt.gz");
//        this.mkSNPsummary_step2("/Users/Aoyue/Documents/chr001_vmap2.1_AnnoDB_10000lines.txt", "/Users/Aoyue/Documents/chr001.wheat.ancestralAllele_1000000.txt", "/Users/Aoyue/Documents/out/chr001_vmap2.1_.txt.gz");
//    this.getCDSannotation("/Users/Aoyue/Documents/test", "/Users/Aoyue/Documents/out");
//        this.classifySNP_byPop("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/102_cdsAnnoDB", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/001_ori");
//        this.changeChrPos();
//        this.mergebySub();
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/001_total", "100", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/002_abd", "100", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/003_ab", "100", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/daf/004_d", "100", "1");
        //10 bins
//        this.mkBarplotofDAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/003_daf/001_mkBarplotofDAF", "10", "1");
//        this.mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/004_maf/001_mkBarplotofMaf","10","0.5");
        //20 bins
//        this.mkBarplotofDAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/005_daf_20bins", "20", "1");
//        this.changeChrPos();
//this.mergebySub();
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/001_genicSNPAnnotation";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/002_changeChrPos";
        /**
         * 新的一批数据，和老师一起写程序得出的结果
         */
//        new CountSites().mergefile1and2_chr1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/002_changeChrPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/003_mergebySub");
//        new CountSites().changechrPosonTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/001_genicSNPAnnotation","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/002_changeChrPos");
//    this.mkBarplot("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/001_test", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/002_mkBarplot", "20", "1.00001");
//    new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/004_merge/chr_A.SNP_anno.txt.gz", "A_SNP_anno.changeChrPos.txt.gz");
//        new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/004_merge/chr_B.SNP_anno.txt.gz", "B_SNP_anno.changeChrPos.txt.gz");
//    new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/003_mergebySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/004_merge/chr_D.SNP_anno.txt.gz", "D_SNP_anno.changeChrPos.txt.gz");
//    new CountSites().mergefile1and2_chr1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/102_cdsAnnoDB","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/001_merge");
//    new CountSites().changechrPosonTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/001_merge","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/002_changeChrPos");
//    new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/002_changeChrPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/003_mergeSub/chr_A.SNP_anno.txt.gz","A_vmap2.1_AnnoDB_addDAF_addSIFT_CDSregion.changeChrPos.txt.gz");
//        new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/002_changeChrPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/003_mergeSub/chr_B.SNP_anno.txt.gz","B_vmap2.1_AnnoDB_addDAF_addSIFT_CDSregion.changeChrPos.txt.gz");
//    new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/002_changeChrPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/105_from102/003_mergeSub/chr_D.SNP_anno.txt.gz","D_vmap2.1_AnnoDB_addDAF_addSIFT_CDSregion.changeChrPos.txt.gz");
        /**
         * l老师新过滤的gene，生成的数据，这里基因数增加
         *
         */
//    new CountSites().mergefile1and2_chr1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/006_genicSNPAnnotation_lowerFilterGene","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/007_changeChrPos");
//    new CountSites().changechrPosonTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/007_changeChrPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/008_mergebySub");
//    new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/008_mergebySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/009_merge/chr_A.SNP_anno.txt.gz","A_SNP_anno.changeChrPos.txt.gz");
//       new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/008_mergebySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/009_merge/chr_B.SNP_anno.txt.gz","B_SNP_anno.changeChrPos.txt.gz");
//        new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/008_mergebySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/009_merge/chr_D.SNP_anno.txt.gz","D_SNP_anno.changeChrPos.txt.gz");
        /**
         * 最新的一批数据，添加了gerp和phyloP的结果
         *
         */
//    new AoMath().filterValue("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/010_genicSNPAnnotation_addGERPandPhyloP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/011_filterGERPandPhylop");
//    new CountSites().mergefile1and2_chr1and2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/010_genicSNPAnnotation_addGERPandPhyloP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/012_changeChrPos");
//        new CountSites().changechrPosonTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/012_changeChrPos","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/013_mergebySub");
//        new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/013_mergebySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_A.SNP_anno.txt.gz","A_SNP_anno.changeChrPos.txt.gz");
//        new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/013_mergebySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_B.SNP_anno.txt.gz","B_SNP_anno.changeChrPos.txt.gz");
//        new CountSites().mergeTxtbysuffix("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/013_mergebySub","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_D.SNP_anno.txt.gz","D_SNP_anno.changeChrPos.txt.gz");
//        this.addSIFTGroup("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/015_addSIFTgroup");
//        new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/015_addSIFTgroupbasedGerpPhyloP","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/015_merge015/chrAll_SNP_anno_group_basedGerpPhyloP.txt.gz");
//        this.addGroup_basedSIFTvalue();
//        new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/016_addSIFTgroupbasedSIFTvalue","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/017_merge016/chrAll_SNP_anno_group_basedSIFTvalue.txt.gz");
//        this.countDeleteriousSNPs_basedCHR();
//        this.countDeleteriousSNPs_basedSubgenome();
//        new AoMath().countValue("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge014");
        /**
         * 计算CDS区域，所有变异的类型个数
         *
         */

//        new CountSites().mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge014/chrAll_anno_original.txt.gz");
//        this.countDeleteriousSNPs_basedSubgenome();

        /**
         * 检查原文件中有没有 同义突变是值小于0.05的   有！！！
         */
//        new AoMath().countValue("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/009_output/000_xls");
//        this.addSIFTGroup_basedGERP();

        /**
         * 计算CDS区域，所有变异的类型个数,这里分组是只有GERP大于1
         *
         */
//        this.countDeleteriousSNPs_basedSubgenome();

        /**
         * 计算原文件中，个体VCF的delterious占 DAF common和rare 的比例
         */

//        this.calIndivDeleteriousProportion();

    }




    public void calIndivDeleteriousProportion(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/022_indivDeleteriousProportion/001_indivcf/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/022_indivDeleteriousProportion/002_proportion/001_hexaploid_deleteriousProportion.txt";
        String annoS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/013_mergebySub/chr1A_SNP_anno.changeChrPos.txt.gz";
        new AoFile().readheader(annoS);
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        try {
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
//            String header = "Individual\tTotalDeleteriousSNPs\tCommonModerate\tCommonLarge\tCommonExtreme\tRareModerate\tRareLarge\tRareExtreme";
            String header = "Indi\tProportion\tEffect";
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();

                BufferedReader br = null;
                if (annoS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(annoS);
                } else if (annoS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(annoS);
                }

                TIntArrayList poslist = new AoFile().getNumListfromVCF(infileS); //引用别的方法
                TDoubleArrayList CommonModerate = new TDoubleArrayList();
                TDoubleArrayList CommonLarge = new TDoubleArrayList();
                TDoubleArrayList CommonExtreme = new TDoubleArrayList();
                TDoubleArrayList RareModerate = new TDoubleArrayList();
                TDoubleArrayList RareLarge = new TDoubleArrayList();
                TDoubleArrayList RareExtreme = new TDoubleArrayList();
                int cntall = 0;
                int cntCommonModerate = 0;
                int cntCommonLarge = 0;
                int cntCommonExtreme = 0;
                int cntRareModerate = 0;
                int cntRareLarge = 0;
                int cntRareExtreme = 0;

                double sift = Double.NaN;
                double gerp = Double.NaN;
                double phylop = Double.NaN;
                double daf = Double.NaN;
                String temp = br.readLine();
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String type = l.get(11);
                    String siftscore = l.get(12);
                    String dafABD = l.get(15);
                    String gerpscore = l.get(17);
                    String phylopscore = l.get(18);
                    int index = poslist.binarySearch(pos);
                    if (index < 0) continue;
                    if (type.equals("SYNONYMOUS"))continue;
                    if (type.equals("NONSYNONYMOUS")) {
                        if (!siftscore.startsWith("N")) {
                            sift = Double.parseDouble(siftscore);
                            if (sift < 0.05) {
                                if (!gerpscore.startsWith("N")) { //均有值存在
                                    gerp = Double.parseDouble(gerpscore);
                                    if (gerp > 1 && gerp < 2.5) {
                                        if (!dafABD.startsWith("N")){
                                            daf = Double.parseDouble(dafABD);
                                            if (daf >= 0.05){//common allele
                                                cntCommonModerate++;
                                                cntall++;
                                            }
                                            else{
                                                cntRareModerate++;
                                                cntall++;
                                            }
                                        }
                                    }
                                    if(gerp >= 2.5 && gerp <3.6){
                                        if (!dafABD.startsWith("N")){
                                            daf = Double.parseDouble(dafABD);
                                            if (daf >= 0.05){//common allele
                                                cntCommonLarge++;
                                                cntall++;
                                            }
                                            else{
                                                cntRareLarge++;
                                                cntall++;
                                            }
                                        }

                                    }
                                    if(gerp >= 3.6){ //max 3.5
                                        if (!dafABD.startsWith("N")){
                                            daf = Double.parseDouble(dafABD);
                                            if (daf >= 0.05){//common allele
                                                cntCommonExtreme++;
                                                cntall++;
                                            }
                                            else{
                                                cntRareExtreme++;
                                                cntall++;
                                            }
                                        }
                                    }
                                } //有值存在
                            }
                        }
                    }
                }
                br.close();
                int all = cntCommonModerate + cntCommonLarge + cntCommonExtreme + cntRareModerate + cntRareLarge + cntRareExtreme;
                double r1 = (double)cntCommonModerate/(double) cntall;
                double r2 = (double)cntCommonLarge/(double) cntall;
                double r3 = (double)cntCommonExtreme/(double) cntall;
                double r4 = (double)cntRareModerate/(double) cntall;
                double r5 = (double)cntRareLarge/(double) cntall;
                double r6 = (double)cntRareExtreme/(double) cntall;

                bw.write(fs[i].getName().split("_")[2] + "\t" + String.format("%.3f",r1) + "\tCommonModerate\n");
                bw.write(fs[i].getName().split("_")[2] + "\t" + String.format("%.3f",r2) + "\tCommonLarge\n");
//                bw.write(fs[i].getName().split("_")[2] + "\t" + String.format("%.3f",r3) + "\tCommonExtreme\n");
                bw.write(fs[i].getName().split("_")[2] + "\t" + String.format("%.3f",r4) + "\tRareModerate\n");
                bw.write(fs[i].getName().split("_")[2] + "\t" + String.format("%.3f",r5) + "\tRareLarge\n");
//                bw.write(fs[i].getName().split("_")[2] + "\t" + String.format("%.3f",r6) + "\tRareExtreme\n");

//                bw.write(fs[i].getName().split("_")[2] + "\t" + cntall + "\t" + String.format("%.3f",r1)+ "\t" + String.format("%.3f",r2)+ "\t" + String.format("%.3f",r3)
//                        + "\t" + String.format("%.3f",r4)+ "\t" + String.format("%.3f",r5)+ "\t" + String.format("%.3f",r6));
//                bw.newLine();
                System.out.println(cntall + "   " + all);

                System.out.println(fs[i].getAbsolutePath() + " is completed at " + outfileS);
                System.out.println(fs[i].getName().split("_")[2] + " have " + cntall + " deleterious SNPs");
            }
            bw.flush();
            bw.close();

        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }



    /**
     * 根据GERP值和PhyloP的值，还有SIFT值判断，再文件最后再添加一列分组信息 Synonymous
     * Nonsynonymous_tolerent Deleterious
     *
     */
    public void addSIFTGroup_basedGERP() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge014";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/021_addSIFTgroupbasedGERP";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_addSIFTgroup_basedGERP.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt.gz", "_addSIFTgroup_basedGERP.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = null;
                if (outfileS.endsWith(".txt")) {
                    bw = IOUtils.getTextWriter(outfileS);
                } else if (outfileS.endsWith(".txt.gz")) {
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }

                double sift = Double.NaN;
                double gerp = Double.NaN;
                double phylop = Double.NaN;
                String temp = null;
                String header = br.readLine();
                bw.write(header + "\tgroup_GERPvalue");
                bw.newLine();
                List<String> l = new ArrayList();
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    //0ID	1Chr	2Pos	3Ref	4Alt	5Major	6Minor	7Maf	8AAF_ABD	9AAF_AB	10Transcript	11Region	12Variant_type	13SIFT_score	14Ancestral	15DAF	16DAF_ABD	17DAF_AB	18Gerp	19PhyloP
                    l = PStringUtils.fastSplit(temp);
                    String type = l.get(11);
                    String siftscore = l.get(12);
                    String gerpscore = l.get(17);
                    String phylopscore = l.get(18);
                    //如果变异类型是同义突变，那么就不用做任何判断；直接加上分组 Synonymous 并写入
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值小于0.05，gerp和phylop存在，且gerp大于1，且phylop大于0.5；那么加上分组 Deleterious 并写入； gerp 值和 phylop值不满足条件的，那么就不进行分组
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值大于0.05，那么加上分组 Nonsynonymous_tolerent 并写入
                    //如果变异类型是非同义突变，SIFT值不存在，那么不分组 不写入
                    if (type.equals("SYNONYMOUS")) {
                        sb.append(temp).append("\tSynonymous");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (type.equals("NONSYNONYMOUS")) {
                        if (!siftscore.startsWith("N")) {
                            sift = Double.parseDouble(siftscore);
                            if (sift < 0.05) {
                                if (!gerpscore.startsWith("N")) { //均有值存在
                                    gerp = Double.parseDouble(gerpscore);
                                    if (gerp > 1) {
                                        sb.append(temp).append("\tDeleterious");
                                        bw.write(sb.toString());
                                        bw.newLine();
                                    }
                                }
                            } else { //sift值大于等于0.05
                                sb.append(temp).append("\tNonsynonymous_tolerant");
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                        }
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }




    /**
     *
     * 根据库文件生成的表格，进行每条染色体，每个类型的计数，输出表格
     */
    public void countDeleteriousSNPs_basedSubgenome() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/015_merge015";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/017_merge016";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/019_countCase";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge014";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/019_countCase";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/021_addSIFTgroupbasedGERP";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/019_countCase";
//        new AoFile().readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge014/chrAll_anno_original.txt.gz");
        new AoFile().readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/021_addSIFTgroupbasedGERP/chrAll_anno_original_addSIFTgroup_basedGERP.txt.gz");
//如何将1A 1B 1D 2A 2B 2D 建立subgenome set联系起来？？？
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 1; i < 8; i++) {
            hm.put(String.valueOf(i) + "A", "Asub");
            hm.put(String.valueOf(i) + "B", "Bsub");
            hm.put(String.valueOf(i) + "D", "Dsub");
        }
        List<String> subgenomeSet = new ArrayList<String>(new HashSet<String>(hm.values()));
        Collections.sort(subgenomeSet);
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_countCase_bySub_basedonlyGERP.txt").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                RowTable<String> t = new RowTable<>(infileS);
                List<String> chr = t.getColumn(0);
                List<String> group = t.getColumn(19);
                List<String> chrSet = new ArrayList<>(new HashSet<>(chr));
                List<String> groupSet = new ArrayList<>(new HashSet<String>(group));
                Collections.sort(groupSet);
                Collections.sort(chrSet);
                int[][] ss = new int[subgenomeSet.size()][groupSet.size()];


                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    String chrp = l.get(0);
                    String sub = hm.get(chrp);
                    String typep = l.get(19);
                    int index1 = Collections.binarySearch(subgenomeSet, sub);
                    int index2 = Collections.binarySearch(groupSet, typep);
                    ss[index1][index2]++;
                }

                //先写表头
                bw.write("Chr");
                for (int i = 0; i < groupSet.size(); i++) {
                    String type = groupSet.get(i);
                    bw.write("\t" + type);
                }
                bw.newLine();
                //再写统计的数目
                for (int i = 0; i < ss.length; i++) {
                    bw.write(subgenomeSet.get(i));
                    for (int j = 0; j < ss[0].length; j++) {
                        String count = Integer.toString(ss[i][j]);
                        bw.write("\t" + count);
                    }
                    bw.newLine();
                }

                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
    }


    /**
     *
     * 根据库文件生成的表格，进行每条染色体，每个类型的计数，输出表格
     */
    public void countDeleteriousSNPs_basedCHR() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/015_merge015";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/019_countCase";
        
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_countCase_byChr.txt").getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_countCase_byChr.txt").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);

                RowTable<String> t = new RowTable<>(infileS);
                List<String> chr = t.getColumn(0);
                List<String> group = t.getColumn(19);
                List<String> chrSet = new ArrayList<>(new HashSet<>(chr));
                List<String> groupSet = new ArrayList<>(new HashSet<String>(group));
                Collections.sort(groupSet);
                Collections.sort(chrSet);
                int[][] ss = new int[chrSet.size()][groupSet.size()];

                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    String chrp = l.get(0);
                    String typep = l.get(19);
                    int index1 = Collections.binarySearch(chrSet, chrp);
                    int index2 = Collections.binarySearch(groupSet, typep);
                    ss[index1][index2]++;
                }

                //先写表头
                bw.write("Chr");
                for (int i = 0; i < groupSet.size(); i++) {
                    String type = groupSet.get(i);
                    bw.write("\t" + type);
                }
                bw.newLine();
                //再写统计的数目
                for (int i = 0; i < ss.length; i++) {
                    bw.write(chrSet.get(i));
                    for (int j = 0; j < ss[0].length; j++) {
                        String count = Integer.toString(ss[i][j]);
                        bw.write("\t" + count);
                    }
                    bw.newLine();
                }

                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }

        });

    }

    /**
     * 根据GERP值和PhyloP的值，还有SIFT值判断，再文件最后再添加一列分组信息 Synonymous
     * Nonsynonymous_tolerent Deleterious
     *
     */
    public void addGroup_basedSIFTvalue() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/016_addSIFTgroupbasedSIFTvalue";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_group_basedSIFTvalue.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt.gz", "_group_basedSIFTvalue.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = null;
                if (outfileS.endsWith(".txt")) {
                    bw = IOUtils.getTextWriter(outfileS);
                } else if (outfileS.endsWith(".txt.gz")) {
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }

                double sift = Double.NaN;
                String temp = null;
                String header = br.readLine();
                bw.write(header + "\tgroup_SIFTvalue");
                bw.newLine();
                List<String> l = new ArrayList();
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    //0ID	1Chr	2Pos	3Ref	4Alt	5Major	6Minor	7Maf	8AAF_ABD	9AAF_AB	10Transcript	11Region	12Variant_type	13SIFT_score	14Ancestral	15DAF	16DAF_ABD	17DAF_AB	18Gerp	19PhyloP
                    l = PStringUtils.fastSplit(temp);
                    String type = l.get(11);
                    String siftscore = l.get(12);
                    //如果变异类型是同义突变，那么就不用做任何判断；直接加上分组 Synonymous 并写入
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值小于0.5，gerp和phylop存在，且gerp大于1，且phylop大于0.5；那么加上分组 Deleterious 并写入； gerp 值和 phylop值不满足条件的，那么就不进行分组
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值大于0.5，那么加上分组 Nonsynonymous_tolerent 并写入
                    //如果变异类型是非同义突变，SIFT值不存在，那么不分组 不写入
                    if (type.equals("SYNONYMOUS")) {
                        sb.append(temp).append("\tSynonymous");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (type.equals("NONSYNONYMOUS")) {
                        if (!siftscore.startsWith("N")) {
                            sift = Double.parseDouble(siftscore);
                            if (sift < 0.05) {
                                sb.append(temp).append("\tDeleterious");
                                bw.write(sb.toString());
                                bw.newLine();
                            } else { //sift值大于等于0.05
                                sb.append(temp).append("\tNonsynonymous_tolerent");
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                        }
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    /**
     * 根据GERP值和PhyloP的值，还有SIFT值判断，再文件最后再添加一列分组信息 Synonymous
     * Nonsynonymous_tolerent Deleterious
     *
     */
    public void addSIFTGroup(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_addSIFTgroup.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt.gz", "_addSIFTgroup.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = null;
                if (outfileS.endsWith(".txt")) {
                    bw = IOUtils.getTextWriter(outfileS);
                } else if (outfileS.endsWith(".txt.gz")) {
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }

                double sift = Double.NaN;
                double gerp = Double.NaN;
                double phylop = Double.NaN;
                String temp = null;
                String header = br.readLine();
                bw.write(header + "\tSIFTgroup");
                bw.newLine();
                List<String> l = new ArrayList();
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    //0ID	1Chr	2Pos	3Ref	4Alt	5Major	6Minor	7Maf	8AAF_ABD	9AAF_AB	10Transcript	11Region	12Variant_type	13SIFT_score	14Ancestral	15DAF	16DAF_ABD	17DAF_AB	18Gerp	19PhyloP
                    l = PStringUtils.fastSplit(temp);
                    String type = l.get(11);
                    String siftscore = l.get(12);
                    String gerpscore = l.get(17);
                    String phylopscore = l.get(18);
                    //如果变异类型是同义突变，那么就不用做任何判断；直接加上分组 Synonymous 并写入
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值小于0.05，gerp和phylop存在，且gerp大于1，且phylop大于0.5；那么加上分组 Deleterious 并写入； gerp 值和 phylop值不满足条件的，那么就不进行分组
                    //如果变异类型是非同义突变，且SIFT值存在，且SIFT值大于0.05，那么加上分组 Nonsynonymous_tolerent 并写入
                    //如果变异类型是非同义突变，SIFT值不存在，那么不分组 不写入
                    if (type.equals("SYNONYMOUS")) {
                        sb.append(temp).append("\tSynonymous");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (type.equals("NONSYNONYMOUS")) {
                        if (!siftscore.startsWith("N")) {
                            sift = Double.parseDouble(siftscore);
                            if (sift < 0.05) {
                                if (!gerpscore.startsWith("N") && (!phylopscore.startsWith("N"))) { //均有值存在
                                    gerp = Double.parseDouble(gerpscore);
                                    phylop = Double.parseDouble(phylopscore);
                                    if (gerp > 1 && (phylop > 0.5)) {
                                        sb.append(temp).append("\tDeleterious");
                                        bw.write(sb.toString());
                                        bw.newLine();
                                    }
                                }
                            } else { //sift值大于等于0.05
                                sb.append(temp).append("\tNonsynonymous_tolerent");
                                bw.write(sb.toString());
                                bw.newLine();
                            }
                        }
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    /**
     * 目的：想通过annotation库文件直接获得作图的表格，免去中间手动建立EXCEL和合并结果
     *
     * @param infileDirS
     * @param outfileDirS
     * @param binNum the number of bins that would be divided
     */
    public void mkBarplot(String infileDirS, String outfileDirS, String binNum, String max) {
        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble(max);
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt", bins + "bins" + ".Table.txt")).getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt.gz", bins + "bins" + ".Table.txt")).getAbsolutePath();
            }
            //先建立bound数组
            double[] bound = new double[bins];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double) length / bins * i;
            }

            //开始计算daf
            double[] daf1 = new double[bins];
            TDoubleArrayList dafList1 = new TDoubleArrayList();
            RowTable<String> t = new RowTable<>(infileS);
            int count = t.getRowNumber();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 15).startsWith("N")) { //DAF值所在的那一列 DAF_ABD
                    continue;
                }
                double value = t.getCellAsDouble(i, 15); //DAF值所在的那一列 DAF_ABD
                dafList1.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf1[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf1.length; i++) {
                daf1[i] = daf1[i] / dafList1.size();
            }
            System.out.println(dafList1.size() + "  size");

            //开始计算daf_ABD
            double[] daf2 = new double[bins];
            TDoubleArrayList dafList2 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 16).startsWith("N")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 16); //DAF值所在的那一列
                dafList2.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf2[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf2.length; i++) {
                daf2[i] = daf2[i] / dafList2.size();
            }
            System.out.println(dafList2.size() + "  size");

            //开始计算daf_AB
            double[] daf3 = new double[bins];
            TDoubleArrayList dafList3 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 17).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 17); //DAF值所在的那一列
                dafList3.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf3[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf3.length; i++) {
                daf3[i] = daf3[i] / dafList3.size();
            }
            System.out.println(dafList3.size() + "  size");

            //开始写出文件
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Daf\tDensity_Total\tDensity_ABD\tDensity_AB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2)).append("\t").append(String.format("%.4f", daf1[i])).append("\t").append(String.format("%.4f", daf2[i])).append("\t").append(String.format("%.4f", daf3[i]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName() + " is completed at " + outfileS);
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
     * @param binNum the number of bins that would be divided
     */
    public void mkBarplotofMAF(String infileDirS, String outfileDirS, String binNum, String max) {
        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble(max);
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt", bins + "bins" + ".Table.txt")).getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt.gz", bins + "bins" + ".Table.txt")).getAbsolutePath();
            }
            //先建立bound数组
            double[] bound = new double[bins];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double) length / bins * i;
            }

            //开始计算daf
            double[] daf1 = new double[bins];
            TDoubleArrayList dafList1 = new TDoubleArrayList();
            RowTable<String> t = new RowTable<>(infileS);
            int count = t.getRowNumber();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 2).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 2); //DAF值所在的那一列
                dafList1.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf1[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf1.length; i++) {
                daf1[i] = daf1[i] / dafList1.size();
            }
            System.out.println(dafList1.size() + "  size");

            //开始计算daf_ABD
            double[] daf2 = new double[bins];
            TDoubleArrayList dafList2 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 3).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 3); //DAF值所在的那一列
                dafList2.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf2[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf2.length; i++) {
                daf2[i] = daf2[i] / dafList2.size();
            }
            System.out.println(dafList2.size() + "  size");

            //开始计算daf_AB
            double[] daf3 = new double[bins];
            TDoubleArrayList dafList3 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 4).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 4); //DAF值所在的那一列
                dafList3.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf3[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf3.length; i++) {
                daf3[i] = daf3[i] / dafList3.size();
            }
            System.out.println(dafList3.size() + "  size");

            //开始写出文件
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Maf\tDensity_Total\tDensity_ABD\tDensity_AB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2)).append("\t").append(String.format("%.4f", daf1[i])).append("\t").append(String.format("%.4f", daf2[i])).append("\t").append(String.format("%.4f", daf3[i]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName() + " is completed at " + outfileS);
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
     * @param binNum the number of bins that would be divided
     */
    public void mkBarplotofDAF(String infileDirS, String outfileDirS, String binNum, String max) {
        int bins = Integer.parseInt(binNum);
        double length = Double.parseDouble(max);
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt", bins + "bins" + ".Table.txt")).getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("txt.gz", bins + "bins" + ".Table.txt")).getAbsolutePath();
            }
            //先建立bound数组
            double[] bound = new double[bins];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double) length / bins * i;
            }

            //开始计算daf
            double[] daf1 = new double[bins];
            TDoubleArrayList dafList1 = new TDoubleArrayList();
            RowTable<String> t = new RowTable<>(infileS);
            int count = t.getRowNumber();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 5).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 5); //DAF值所在的那一列
                dafList1.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf1[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf1.length; i++) {
                daf1[i] = daf1[i] / dafList1.size();
            }
            System.out.println(dafList1.size() + "  size");

            //开始计算daf_ABD
            double[] daf2 = new double[bins];
            TDoubleArrayList dafList2 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 6).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 6); //DAF值所在的那一列
                dafList2.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf2[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf2.length; i++) {
                daf2[i] = daf2[i] / dafList2.size();
            }
            System.out.println(dafList2.size() + "  size");

            //开始计算daf_AB
            double[] daf3 = new double[bins];
            TDoubleArrayList dafList3 = new TDoubleArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, 7).equals("NA")) { //DAF值所在的那一列
                    continue;
                }
                double value = t.getCellAsDouble(i, 7); //DAF值所在的那一列
                dafList3.add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                daf3[index]++; //值落入第i种变异的第index个区间的个数
            }
            //开始计算每个区间落入点的比例
            for (int i = 0; i < daf3.length; i++) {
                daf3[i] = daf3[i] / dafList3.size();
            }
            System.out.println(dafList3.size() + "  size");

            //开始写出文件
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Daf\tDensity_Total\tDensity_ABD\tDensity_AB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2)).append("\t").append(String.format("%.4f", daf1[i])).append("\t").append(String.format("%.4f", daf2[i])).append("\t").append(String.format("%.4f", daf3[i]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    public void mergebySub() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }
        for (int i = 0; i < fs.length; i++) {
            new CountSites().mergeTxtbysuffix(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath() + "_A.txt.gz", "A.");
            new CountSites().mergeTxtbysuffix(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath() + "_B.txt.gz", "B.");
            new CountSites().mergeTxtbysuffix(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath() + "_D.txt.gz", "D.");
        }

    }

    public void changeChrPos() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/001_ori";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();

        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }

        for (int i = 0; i < fs.length; i++) {
            new CountSites().mergefileandChangeChrPos_chr1and2(fs[i].getAbsolutePath(), new File(outfileDirS, fs[i].getName()).getAbsolutePath());
        }
    }

    /**
     * 将SNP按照同义非同义突变进行分类，并画出DAF分布图
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void classifySNP_byPop(String infileDirS, String outfileDirS) {

        String[] snpClass = {"delSNP", "nonsyTolerantSNP", "synSNP"};
        String[] out = new String[snpClass.length];
        for (int i = 0; i < snpClass.length; i++) {
            new File(outfileDirS, snpClass[i]).mkdirs();
            out[i] = new File(outfileDirS, snpClass[i]).getAbsolutePath();
        }

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            if (fs[i].getName().endsWith(".xlsx")) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String MAF_ABD = null;
                String MAF_AB = null;
                String chrS = f.getName().substring(3, 6);
                boolean ifd = false;
                //根据染色体号进行AB还是D的判断
                String[] db = {"5", "6", "11", "12", "17", "18", "23", "24", "29", "30", "35", "36", "41", "42"};
                Arrays.sort(db);
                if (Arrays.binarySearch(db, chrS) > -1) { //说明是属于D的
                    ifd = true;
                }

                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }

                BufferedWriter[] bw = new BufferedWriter[snpClass.length];
                for (int i = 0; i < bw.length; i++) {
                    outfileS = new File(out[i], "chr" + chrS + "." + snpClass[i] + ".txt.gz").getAbsolutePath();
                    bw[i] = IOUtils.getTextGzipWriter(outfileS);

                    if (ifd == false) {
                        bw[i].write("Chr\tPos\tMaf\tMAF_ABD\tMAF_AB\tDaf\tDaf_ABD\tDaf_AB\tTrans");
                        bw[i].newLine();
                    } else if (ifd == true) {
                        bw[i].write("Chr\tPos\tMaf\tMAF_ABD\tMAF_D\tDaf\tDaf_ABD\tDaf_D\tTrans");
                        bw[i].newLine();
                    }
                }

                String temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
//0Chr	1Pos	2Ref	3Alt	4Major	5Minor	6Maf	7AAF_ABD	8AAF_AB	9Ancestral	10Daf	11Daf_ABD	12Daf_AB	13Variant_type	14SIFT_score	15Transcript
                    List<String> l = PStringUtils.fastSplit(temp);
                    String chr = l.get(0);
                    String pos = l.get(1);
                    String major = l.get(4);
                    String minor = l.get(5);
                    String maf = l.get(6);
                    String AAF_ABD = l.get(7); //注意AAF_ABD中含有0 和1
                    String AAF_AB = l.get(8);
                    String anc = l.get(9);
                    String daf = l.get(10);
                    String Daf_ABD = l.get(11); //注意Daf_ABD中 不含有0和1
                    String Daf_AB = l.get(12);
                    String snpType = l.get(13);
                    String scoreS = l.get(14);
                    String trans = l.get(15);

                    if (AAF_ABD.equals("0.0000") || AAF_ABD.equals("1.0000")) {
                        MAF_ABD = "NA";
                    }
                    if (AAF_AB.equals("0.0000") || AAF_AB.equals("1.0000")) {
                        MAF_AB = "NA";
                    }
                    if (!AAF_ABD.equals("0.0000") && (!AAF_ABD.equals("1.0000"))) { //AAF有分离
                        MAF_ABD = String.valueOf(Math.min(Double.parseDouble(AAF_ABD), 1 - Double.parseDouble(AAF_ABD)));
                    }
                    if (!AAF_AB.equals("0.0000") && (!AAF_AB.equals("1.0000"))) { //AAF有分离
                        MAF_AB = String.valueOf(Math.min(Double.parseDouble(AAF_AB), 1 - Double.parseDouble(AAF_AB)));
                    }

                    //先过滤没有type类型的位点，只保留有类型的位点
                    if (snpType.equals("NA")) {
                        continue;
                    }
                    if (snpType.equals("NONSYNONYMOUS")) { //在类型下进行sift值的判断，
                        if (scoreS.equals("NA")) {
                            continue;
                        }
                        double score = Double.parseDouble(l.get(14));
                        if (score < 0.05) {//说明是有害突变
                            bw[0].write(chr + "\t" + pos + "\t" + maf + "\t" + MAF_ABD + "\t" + MAF_AB + "\t" + daf + "\t" + Daf_ABD + "\t" + Daf_AB + "\t" + trans);
                            bw[0].newLine();
                        } else {//说明是可忍受突变
                            bw[1].write(chr + "\t" + pos + "\t" + maf + "\t" + MAF_ABD + "\t" + MAF_AB + "\t" + daf + "\t" + Daf_ABD + "\t" + Daf_AB + "\t" + trans);
                            bw[1].newLine();
                        }
                    }
                    if (snpType.equals("SYNONYMOUS")) {
                        bw[2].write(chr + "\t" + pos + "\t" + maf + "\t" + MAF_ABD + "\t" + MAF_AB + "\t" + daf + "\t" + Daf_ABD + "\t" + Daf_AB + "\t" + trans);
                        bw[2].newLine();
                    }
                }
                for (int i = 0; i < snpClass.length; i++) {
                    bw[i].flush();
                    bw[i].close();
                }
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileDirS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getCDSannotation(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "chr");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_CDSregion.txt.gz").getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }

            try {
                String header = br.readLine();
                bw.write(header);
                bw.newLine();
                String temp = null;
                List<String> l = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String trans = l.get(15);
                    if (trans.equals("NA")) {
                        continue;
                        /*如果siftascore的值为NA，则无法判断其为有害或是中性突变。我们要筛选即有sift变异类型又有sift值的sites*/
                    }
                    cnt++;
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\t" + String.valueOf(cnt) + "\t" + "trans sites is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * Goal:将ancestral allele添加到数据库中，并计算Daf,Daf_ABD Daf_AB Daf_D
     *
     * @param dbfileS
     * @param ancS
     * @param outfileS
     */
    public void mkSNPsummary_step2(String dbfileS, String ancS, String outfileS) {
        boolean ifd = false;
        double daf = Double.NaN;
        double daf_ABD = Double.NaN;
        double daf_AB = Double.NaN;
        double daf_D = Double.NaN;
        int cntAncNum = 0;
        File f = new File(ancS); //根据ancestral allele 文件，得到染色体号
        int chr = Integer.parseInt(f.getName().substring(3, 6));
        //根据染色体号进行AB还是D的判断
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, chr) > -1) { //说明是属于D的
            ifd = true;
        }
        TIntArrayList snpPosList = new TIntArrayList();
        HashMap<Integer, String> hm = new HashMap<>();

        /*==================================== 建立ancestral allele HashMap =============================================*/
        try { // chr001.wheat.ancestralAllele.txt  chr001_vmap2.1_AnnoDB.txt.gz
            BufferedReader br = null;
            if (f.getName().endsWith(".txt")) {
                br = IOUtils.getTextReader(ancS);
            } else if (f.getName().endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(ancS);
            }
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                cntAncNum++;
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String anc = l.get(3); //
                snpPosList.add(pos);
                hm.put(pos, anc);
            }
            br.close();
            System.out.println(f.getName() + "\tis completed on posList DB with ancestral allele number " + cntAncNum);
        } catch (Exception e) {
            e.printStackTrace();
        }

        int[] snpPos = snpPosList.toArray(new int[snpPosList.size()]);
        Arrays.sort(snpPos);

        try {
            String chrS = PStringUtils.getNDigitNumber(3, chr);
            BufferedReader br = null;
            if (dbfileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(dbfileS);
            } else if (dbfileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(dbfileS);
            }
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }

            String temp = br.readLine(); //read header

            if (ifd == false) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_AB");
                bw.newLine();
            } else if (ifd == true) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_D");
                bw.newLine();
            }

            int cntAnc = 0;
            int cntAncNotMajororMinor = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String major = l.get(4);
                String minor = l.get(5);
                double maf = Double.parseDouble(l.get(6));
                double AAF_ABD = Double.parseDouble(l.get(7));
                double AAF_AB = Double.parseDouble(l.get(8));
                int index = Arrays.binarySearch(snpPos, pos);
                StringBuilder sb = new StringBuilder();
                if (index > -1) { //表明含有anc
                    String ancAllele = hm.get(pos);
                    //如果ancestral allele存在,且等于major，则derived allele等于minor, daf 就等于maf
                    //如果ancestral allele存在,且等于minor，则derived allele等于major, daf 就等于 1-daf1
                    if (ancAllele.equals(minor)) {
                        cntAnc++;
                        daf = 1 - maf;
                        if (AAF_ABD > 0.5) { //说明AAF_ABD是major， 祖先状态是minor的，所有DAF是major
                            daf_ABD = AAF_ABD;
                        } else if (AAF_ABD < 0.5) { //说明AAF_ABD是minor， 祖先状态是minor的，所有DAF是major
                            daf_ABD = 1 - AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = 1 - AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (ancAllele.equals(major)) {
                        cntAnc++;
                        daf = maf;
                        if (AAF_ABD > 0.5) { //说明AAF_ABD是major， 祖先状态是major的，所有DAF是minor
                            daf_ABD = 1 - AAF_ABD;
                        } else if (AAF_ABD < 0.5) {
                            daf_ABD = AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = 1 - AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (!ancAllele.equals(minor) && (!ancAllele.equals(major))) {
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                        cntAncNotMajororMinor++;
                    }

                } else { //表明不含anc
                    sb.append(temp).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            double ratio = (double) cntAncNotMajororMinor / (cntAncNotMajororMinor + cntAnc);
            bw.flush();
            bw.close();
            br.close();
            System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
            System.out.println(new File(dbfileS).getName() + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor. The ratio is " + String.format("%.4f", ratio));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 目的：1.将vmap2的chr pos 提取出来，建立数据库；
     *
     * @param infileS
     * @param outfileS
     */
    public void mkSNPsummary_step1(String infileS, String outfileS) {
        //Chr	Pos	Ref	Alt	Major	Minor	Maf	AAF_ABD	AAF_AB
        boolean ifd = false;
        File f = new File(infileS);
        int CHR = Integer.parseInt(f.getName().substring(3, 6));
        //根据染色体号进行AB还是D的判断
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, CHR) > -1) { //说明是属于D的
            ifd = true;
        }

        try {
            BufferedReader br = null;
            BufferedWriter bw = null; // IOUtils.getTextGzipWriter(outfileS);
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

            if (ifd == false) {
                bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB");
                bw.newLine();
            } else if (ifd == true) {
                bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D");
                bw.newLine();
            }

            String temp = null;
            String te[] = null;
            String major = null;
            String minor = null;
            int biallelicNum = 0;

            String AAF_ABD = null;
            String AAF_AB = null;
            String AAF_D = null;

            while ((temp = br.readLine()) != null) {
                int genoNum = 0;
//                double homNum = 0;
//                double hetNum = 0;
//                double hetRate = 0;
//                double missNum = 0;
//                double missRate = 0;

                double refAlleleGametes = 0;
                double altAlleleGametes = 0;
                double refAF = 0;
                double altAF = 0;
                double maf = 0;
                if (temp.startsWith("#")) {
                    //bw.write(temp);
                    //bw.newLine();
                } else {
                    te = temp.split("\t");
                    String chr = PStringUtils.fastSplit(temp).get(0);
                    String pos = PStringUtils.fastSplit(temp).get(1);
                    String ref = PStringUtils.fastSplit(temp).get(3);
                    String alt = PStringUtils.fastSplit(temp).get(4);

                    if (te[4].length() == 1) { //不含有逗号的情况，即只有一个alt。又开始分，是D 是I 是ATGC 3种情况
                        if (alt.contains("D") || alt.contains("I")) {
                            continue; //只有一个alt且不是indel
                        }
                        AAF_ABD = te[7].split(";")[7].split("=")[1];
                        AAF_AB = te[7].split(";")[8].split("=")[1];
                        AAF_D = te[7].split(";")[8].split("=")[1];

                        biallelicNum++;
                        for (int i = 9; i < te.length; i++) {
                            if (te[i].startsWith(".")) {
//                                missNum++;
                            }
                            if (!te[i].startsWith(".")) {
                                genoNum++; //have the genotype
                                if (te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
//                                    hetNum++; //the number of heterozygous
                                    refAlleleGametes++;
                                    altAlleleGametes++;
                                }
                                if (te[i].startsWith("0/0")) {
//                                    homNum++; //the number of heterozygous
                                    refAlleleGametes++;
                                    refAlleleGametes++;
                                }
                                if (te[i].startsWith("1/1")) {
//                                    homNum++;
                                    altAlleleGametes++;
                                    altAlleleGametes++;
                                }
                            }
                        }
//                        hetRate = hetNum / genoNum;
//                        missRate = missNum / (missNum + genoNum);
                        refAF = refAlleleGametes / (refAlleleGametes + altAlleleGametes);
                        altAF = altAlleleGametes / (refAlleleGametes + altAlleleGametes);;
                        if (refAF > altAF) {
                            major = ref;
                            minor = alt;
                            maf = altAF;
                        } else {
                            maf = refAF;
                            major = alt;
                            minor = ref;
                        }
                        StringBuilder sb = new StringBuilder();
                        //bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB");
                        sb.append(chr).append("\t").append(pos).append("\t").append(ref).append("\t").append(alt).append("\t").
                                append(major).append("\t").append(minor).append("\t").append(String.format("%.4f", maf)).append("\t").
                                append(AAF_ABD).append("\t").append(AAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();

                    }
                }
            }
            System.out.println(infileS + " is completed at " + outfileS);
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 对文件进行分bin，画分布图
     */
    private void mkBarplotOfSNPs() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/delSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/delSNPCount.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/delSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/delSNP/dafSFS.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/nonsyTolerantSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/nonsyTolerantSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/nonsyTolerantSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/nonsyTolerantSNP/dafSFS.txt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/004_merge/synSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/synSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/synSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/005_bin/synSNP/dafSFS.txt";
        /**
         * ******************** D subgenome ****************************
         */
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/delSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/delSNPCount.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/delSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/delSNP/dafSFS.txt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/nonsyTolerantSNP/";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/nonsyTolerantSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/nonsyTolerantSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/nonsyTolerantSNP/dafSFS.txt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/009_merge/synSNP";
//        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/synSNP.txt";
//        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/synSNP/mafSFS.txt";
//        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/006_snpClassify/010_bin/synSNP/dafSFS.txt";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub";
        String countFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/synSNP.txt";
        String mafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/mafSFS.txt";
        String dafDistrubutionFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/004_bin/dafSFS.txt";

        //int sampleSize = 10000;
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                //System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);//将文件删除后，重新将文件列表打印出来，此时，fs不包含隐藏文件。
        /*建立一个边界数组bound，大小是100；bound[i] = 1/100*i;即，均分为100等分！
        建立一个二维数组mafFrequency 和 dafFrequency， 长度为class的种类长，100宽；
        建立一个count 和 dafCount 数组，长度为class分类长；
        建立一个 dafList1 和 dafList 数组，长度为class分类长；
        
        进入for循环，对class文件一一遍历，以读表格的形式进入文件
        Chr	Pos	MinorAllele	MAF	DerivedAllele	DAF
        1	92716	A	0.0077619664	NA	NA
        1	93774	G	8.130081E-4	A	0.9991869919
        1	122123	C	0.0012254902	C	0.0012254902
        做以下几件事情：1，数行数，看每个分类的个数；每个位点一定有maf的值，但不一定有daf的值，因为DA allele不一定存在，若不存在，则DAF值为空。
        2，将maf的值加入mafList;在bound数组里搜索，如果未搜到，则index = - index -2; 否则mafFrequency[i][index]++;
        如果DAF不以N开头，则将daf的值加入 dafList;在bound数组里搜索，如果未搜到，则index = - index -2; 否则dafFrequency[i][index]++; dafCount数组加一
        在此循环内，进入进入另一个for循环j，做统计：
        计算出maf 和daf 在 index 1-100 bin范围内,各个位点所占的百分比。即 index=1, 所有位点count= rownumber， frenquency = index /count；
         */
        int size = 100; //把daf值分成100份，每份有1/100=0.001长度
        double[] bound = new double[size];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double) 1 / size * i;
        }
        double[][] mafFrequency = new double[fs.length][size]; //变异类型的个数
        double[][] dafFrequency = new double[fs.length][size];
        int[] count = new int[fs.length]; //每个变异类型的个数
        int[] dafCount = new int[fs.length];
        TDoubleArrayList[] mafList = new TDoubleArrayList[fs.length]; //每种变异类型的值的集合
        TDoubleArrayList[] dafList = new TDoubleArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            mafList[i] = new TDoubleArrayList();
            dafList[i] = new TDoubleArrayList();
            String infileS = fs[i].getAbsolutePath();
            RowTable<String> t = new RowTable<>(infileS);
            count[i] = t.getRowNumber();
            for (int j = 0; j < t.getRowNumber(); j++) {
                double value = t.getCellAsDouble(j, 2);
                mafList[i].add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                mafFrequency[i][index]++; //值落入第i种变异的第index个区间的个数
                if (!t.getCell(j, 3).startsWith("N")) {
                    value = t.getCellAsDouble(j, 3);
                    dafList[i].add(value);
                    index = Arrays.binarySearch(bound, value);
                    if (index < 0) {
                        index = -index - 2;
                    }
                    dafFrequency[i][index]++;
                    dafCount[i]++;
                }
            }
            for (int j = 0; j < mafFrequency[i].length; j++) {
                mafFrequency[i][j] = mafFrequency[i][j] / count[i];
                dafFrequency[i][j] = dafFrequency[i][j] / dafCount[i]; //因为daf值不是每个都有，有些pos是NA值，所以需要重新计算。
            }
        }
        /*打表格，输入表头，去掉最后一个\t键；表头为4个文件的文件名；
        第二行输入每个分类的conut数；
         */
        try {
            BufferedWriter bw = IOUtils.getTextWriter(countFileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fs.length; i++) {
                sb.append(fs[i].getName().replaceFirst(".changeChrPos.txt", "")).append("\t");
            }
            sb.deleteCharAt(sb.length() - 1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length - 1; i++) {
                bw.write(String.valueOf(count[i]) + "\t");
            }
            bw.write(String.valueOf(count[count.length - 1]));
            bw.newLine();
            bw.flush();
            bw.close();

            bw = IOUtils.getTextWriter(mafDistrubutionFileS); //开始写新的文件MAF
            bw.write("MAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t" + fs[i].getName().replaceFirst(".changeChrPos.txt", "")); //循环写表头
            }
            bw.newLine();
            /*double[][] mafFrequency = new double[fs.length][size] 文件长度为4，size为100*/
            for (int i = 0; i < mafFrequency[0].length; i++) { //i小于第一个文件的长度100，
                bw.write(String.format("%.2f", bound[i]));
                for (int j = 0; j < mafFrequency.length; j++) { //将1-100的频率写出来
                    bw.write("\t" + mafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();

            bw = IOUtils.getTextWriter(dafDistrubutionFileS);//开始写新的文件DAF
            bw.write("DAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t" + fs[i].getName().replaceFirst(".changeChrPos.txt", ""));
            }
            bw.newLine();
            for (int i = 0; i < dafFrequency[0].length; i++) {
                bw.write(String.format("%.2f", bound[i]));
                for (int j = 0; j < dafFrequency.length; j++) {
                    bw.write("\t" + dafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 将SNP按照同义非同义突变进行分类，并画出DAF分布图
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void classifySNPs(String infileDirS, String outfileDirS) {
        String[] snpClass = {"delSNP", "nonsyTolerantSNP", "synSNP"};
        String[] out = new String[snpClass.length];
        for (int i = 0; i < snpClass.length; i++) {
            new File(outfileDirS, snpClass[i]).mkdirs();
            out[i] = new File(outfileDirS, snpClass[i]).getAbsolutePath();
        }

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
                String chrS = f.getName().substring(3, 6);
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }

                BufferedWriter[] bw = new BufferedWriter[snpClass.length];
                for (int i = 0; i < bw.length; i++) {
                    outfileS = new File(out[i], "chr" + chrS + "." + snpClass[i] + ".txt.gz").getAbsolutePath();
                    bw[i] = IOUtils.getTextGzipWriter(outfileS);
                    bw[i].write("Chr\tPos\tMaf\tDaf\tTrans");
                    bw[i].newLine();
                }
                String temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    String chr = l.get(0);
                    String pos = l.get(1);
                    String major = l.get(4);
                    String minor = l.get(5);
                    double maf = Double.parseDouble(l.get(6));
                    String snpType = l.get(9);
                    String scoreS = l.get(10);
                    String trans = l.get(11);
                    String anc = l.get(12);
                    String daf = l.get(13);
                    //先过滤没有type类型的位点，只保留有类型的位点
                    if (snpType.equals("NA")) {
                        continue;
                    }
                    if (snpType.equals("NONSYNONYMOUS")) { //在类型下进行sift值的判断，
                        if (scoreS.equals("NA")) {
                            continue;
                        }
                        double score = Double.parseDouble(l.get(10));
                        if (score < 0.05) {//说明是有害突变
                            bw[0].write(chr + "\t" + pos + "\t" + maf + "\t" + daf + "\t" + trans);
                            bw[0].newLine();
                        } else {//说明是可忍受突变
                            bw[1].write(chr + "\t" + pos + "\t" + maf + "\t" + daf + "\t" + trans);
                            bw[1].newLine();
                        }
                    }
                    if (snpType.equals("SYNONYMOUS")) {
                        bw[2].write(chr + "\t" + pos + "\t" + maf + "\t" + daf + "\t" + trans);
                        bw[2].newLine();
                    }
                }
                for (int i = 0; i < snpClass.length; i++) {
                    bw[i].flush();
                    bw[i].close();
                }
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileDirS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 解析玉米的基因总结分析
     */
    private void summarizeTranscript2() {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
//        String infileDirS = "/data1/home/aoyue/maizeLoad/001_variantSummary/001_hmp321Info_filter";
//        String geneFeatureFileS = "/data1/publicData/maize/gene/Zea_mays.AGPv4.38.pgf";
//        String outfileS = "/data1/home/aoyue/maizeLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
        double gerpCut = 0;
        File[] fs = new File(infileDirS).listFiles();
        int chrNum = fs.length;
        HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
        int[][] snpPos = new int[chrNum][];
        byte[][] snps = new byte[chrNum][];
        byte[][] snpAnc = new byte[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            posGeneMap[i] = new HashMap();
        }
        //下面这一段将posGeneMap建立完整，使每个位点对应哪些基因名字，都装进这个map里 
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        HashMap<String, Integer> geneCDSLengthMap = new HashMap();
        /*将所有基因的名字进行for循环输入到数组genes中，对应于每一个基因，我们通过getTranscriptName得到转录本的名字，通过getCDSList方法得到编码序列的起始位点*/
        List<String> genesList = new ArrayList<>();
        //String[] genes = new String[gf.getGeneNumber()]; //这个是原来的代码，后续需要修改，因为gf文件中包含我们不需要的染色体上的基因
        int cntchr11and12 = 0;
        int cntchr1to10 = 0;

        //*********************************** START1 ***********************************//
        //该段代码的作用是，通过读取每个基因，得到最长转录本的名字，计算该转录本的长度。
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chrIndex = gf.getGeneChromosome(i) - 1;
            /*这个地方是先过滤数据，将定位在11号12号染色体上的基因过滤掉，并且跳出循环*/
            if (chrIndex > 9) {
                cntchr11and12++;
                continue;
            }
            cntchr1to10++; //能够得到1-10号染色体的基因数目
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String geneName = gf.getTranscriptName(i, longTransIndex); //得到最长的转录本的名字
            //genes[i] = geneName;
            genesList.add(geneName);
            List<Range> cdsList = gf.getCDSList(i, longTransIndex);
            /*得到基因的最长转录本的CDSList*/
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
                    } else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                        /*最终将posGeneMap绘图完成*/
                    }
                    cnt++;
                    /*每一个CDS位点相加，最终得到这个cds的长度。*/
                }

                // 最终cnt是一个基因的所有cdslist中，每个cds的每个位点包含的基因数目的总和
            } //该循环是一个基因的所有cds循环
            geneCDSLengthMap.put(geneName, cnt); //
        }
        //*********************************** END1 ***********************************//

        System.out.println(cntchr11and12 + "genes are not used");
        System.out.println(cntchr1to10 + "genes are used");

        String[] genes = genesList.toArray(new String[genesList.size()]);
        Arrays.sort(genes); //genes指所有基因对应最长转录本的名字的组合，是一个数组。
        int[] snpCount = new int[genes.length]; //这一步也很重要，就是想知道每个基因内部有多少个变异位点。
        List<File> hmpList = Arrays.asList(fs); //fs指的是根据VCF位点建立的注释数据库，该数据库包含很多计算的信息，有SIFT值GERP值。
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", "")) - 1;
            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            TByteArrayList snpAncList = new TByteArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                /*这一段主要是将有基因的位点列出来，及将转录组的位点列出来*/
                while ((temp = br.readLine()) != null) {
                    cnt++; //检测程序运行的情况
                    if (cnt % 1000000 == 0) {
                        System.out.println("Hmp\tchr" + String.valueOf(chrIndex + 1) + "\t" + String.valueOf(cnt) + " ###hmpInfo Process");
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).contains("<") || l.get(3).contains(",")) {
                        continue; //过滤含有2个alt和含有indel的位点
                    }
                    int pos = Integer.valueOf(l.get(1));
                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) {
                        continue; //说明该变异位点不在基因区
                    }
                    for (int i = 0; i < geneNameList.size(); i++) {
                        int index = Arrays.binarySearch(genes, geneNameList.get(i));//在基因库的第i个位置，该基因数加一
                        snpCount[index]++; //第i个位置的基因含有的snp数目；
                    }
                    /*将该位置加入snpPosList，表明该位置是有基因的，没有基因的位置就不必加到snpPosList中去*/
                    snpPosList.add(pos);
                    snpList.add(l.get(3).getBytes()[0]); //返回该位点的alt碱基的二进制码
                    snpAncList.add(l.get(4).getBytes()[0]);  //返回该位点的ancestral碱基的二进制码
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            snpPos[chrIndex] = snpPosList.toArray(); //第n条染色体的在基因区间的snp所对应的pos的集合；
            snps[chrIndex] = snpList.toArray(); //第n条染色体的在区间内的snp所对应的alt的集合
            snpAnc[chrIndex] = snpAncList.toArray(); //第n条染色体的在区间内的snp所对应的ancestral allele的集合
        });

        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length]; //非同义突变的数目
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];  //非同义突变，但是sift值是NA的数目

        int[] b73SynCount = new int[genes.length];
        int[] b73NonCount = new int[genes.length];
        int[] b73DelCount = new int[genes.length];
        int[] b73DelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];

        TIntArrayList[] delPosList = new TIntArrayList[chrNum]; //有害变异的位点集合
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", "")) - 1;
            delPosList[chrIndex] = new TIntArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt % 1000000 == 0) {
                        System.out.println("Sift\tchr" + String.valueOf(chrIndex + 1) + "\t" + String.valueOf(cnt) + " ### SIFT Process");
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(16).startsWith("NA")) {
                        continue; //Variant_type	SIFT_score	Transcript  16 17 18
                    }
                    if (l.get(18).startsWith("NA")) {
                        continue; //没有变异类型和转录本的位点，都过滤掉。
                    }
                    String gene = l.get(18);
                    int geneIndex = Arrays.binarySearch(genes, gene); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) {
                        continue;
                    }
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos); // 
                    if (index < 0) {
                        continue;
                    }
                    if (snps[chrIndex][index] != l.get(3).getBytes()[0]) {
                        continue; //再次验证，如果该位点的alt和数据库中的alt不一致，则过滤掉。
                    }
                    byte ref = l.get(2).getBytes()[0];
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (snpAnc[chrIndex][index] == snps[chrIndex][index]) { //如果ancestral allele 等于alt的话，derived allele就等于1
                        derivedState = 1; //mean b73 carries derived allele
                    } else if (snpAnc[chrIndex][index] == ref) { //如果ancestral allele 等于ref的话，derived allele就等于0
                        derivedState = 0;
                    }

                    if (derivedState == -1) {
                        noAncCount[geneIndex]++; //如果ancestral allele 不存在的话，derived allele就等于-1
                    }
                    String type = null;
                    if (l.get(16).equals("NA")) {

                    } else {
                        if (l.get(16).equals("SYNONYMOUS")) { //如果type等于syn，那么 该位点所属的基因的syn属性就加一
                            type = "Syn";
                            synCount[geneIndex]++;
                            if (derivedState == 1) {  //如果derived allele就等于1，说明ancestral allele 等于alt， derived allele 等于ref； 如何计算daf,判断da是major还是minor，如果da是major那么daf=1-daf1，如何判断major allele和minor allele？ 如果ref allele frequency > alt allele frequence,那么major是ref; 反之亦然；
                                b73SynCount[geneIndex]++; //如果参考基因组是 derived allele,那么就加一，为什么？
                            }
                        } else {
                            type = "Non";
                            nonCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73NonCount[geneIndex]++;
                            }
                            if (l.get(17).startsWith("NA")) { //index 17列是sift的值，NON-SYNONYMOUS 存在的情况下，sift可能有也可能没有。
                                naCount[geneIndex]++; //是nonsynonymous类型但是没有sift值的个数
                            } else {
                                if (Double.valueOf(l.get(17)) < 0.05) {
                                    delCount[geneIndex]++;
                                    delPosList[chrIndex].add(pos);
                                    if (derivedState == 1) {
                                        b73DelCount[geneIndex]++;
                                    }
                                }
                            }
                        }
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        /**
         *
         */
        int[][] delPos = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delPos[i] = delPosList[i].toArray();
            Arrays.sort(delPos[i]);
        }
        int[] gerpAlignCount = new int[genes.length];
        int[] snpGerpAlignCount = new int[genes.length];
        double[] gerpTree = new double[genes.length];
        double[] gerpScore = new double[genes.length];
        double[] snpGerpTree = new double[genes.length];
        double[] snpGerpScore = new double[genes.length];

        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", "")) - 1;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) { //gerp文件没有表头
                    cnt++;
                    if (cnt % 1000000 == 0) {
                        System.out.println("Gerp\tchr" + String.valueOf(chrIndex + 1) + "\t" + String.valueOf(cnt) + " ### Gerp Process");
                    }
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(1));
                    if (l.get(14).equals("NA")) {
                        continue;
                    }
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos); //根据pos信息，得到该pos对应的gene name的集合
                    if (geneNameList == null) {
                        continue;
                    }
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        if (gene == null) {
                            continue;
                        }
                        int geneIndex = Arrays.binarySearch(genes, gene);

                        double treeValue = Double.valueOf(l.get(14));
                        double scoreValue = Double.valueOf(l.get(15));
                        if (treeValue == 0) {
                            continue; //过滤枝长是0的数目
                        }
                        gerpAlignCount[geneIndex]++; //如果枝长不是0，说明该位点存在保守不保守
                        gerpTree[geneIndex] += treeValue;
                        gerpScore[geneIndex] += scoreValue; //第i个基因的gerpscore的总和是多少
                        int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                        if (index < 0) {
                            continue;
                        }

                        byte derivedState = 0; //mean ancestral allele is not defined or ancestral allele is alt
                        if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                            derivedState = 1; //mean b73 carries derived allele
                        }
                        snpGerpAlignCount[geneIndex]++;
                        snpGerpTree[geneIndex] += treeValue;
                        snpGerpScore[geneIndex] += scoreValue;
                        index = Arrays.binarySearch(delPos[chrIndex], pos);
                        if (index < 0) {
                            continue;
                        }
                        if (scoreValue <= gerpCut) {
                            continue;
                        }
                        delHGCount[geneIndex]++;
                        if (derivedState == 1) {
                            b73DelHGCount[geneIndex]++;
                        }
                    }

                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        gf.sortGeneByStartPosition();
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpTreeLength\tMeanGerpScore\tSNPMeanGerpTreeLength\tSNPMeanGerpScore";
            header = header + "\tNumAmbigousAnc\tB73NumberOfSyn\tB73PercentageSyn\tB73NumberOfNon\tB73PercentageNon\tB73NumberOfDeleterious\tB73PercentageDeleterious\tB73NumberOfHGDeleterious\tB73PercentageHGDeleterious";
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < genes.length; i++) {
                //if(gf.getGeneChromosome(i) > 10) continue;
                StringBuilder sb = new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append((double) snpCount[i] / cdsLength).append("\t");
                int ifSiftAligned = 1;
                if (naCount[i] == nonCount[i]) {
                    ifSiftAligned = 0; //非同义突变，但是sift值是NA的数目 等于 非同义突变的数目，那么说明在这个基因内部没有sift突变
                }
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append((double) synCount[i] / cdsLength).append("\t");
                sb.append(nonCount[i]).append("\t").append((double) nonCount[i] / cdsLength).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) {
                    ratio = Double.NaN;
                } else {
                    ratio = (double) nonCount[i] / synCount[i];
                }
                sb.append(ratio).append("\t").append(delCount[i]).append("\t").append((double) delCount[i] / cdsLength).append("\t");
                sb.append(delHGCount[i]).append("\t").append((double) delHGCount[i] / cdsLength).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) {
                    ifGerpAligned = 0;
                }
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append((double) gerpAlignCount[i] / cdsLength).append("\t");

                if (gerpAlignCount[i] == 0) {
                    sb.append(Double.NaN).append("\t").append(Double.NaN).append("\t");
                } else {
                    sb.append((double) gerpTree[i] / gerpAlignCount[i]).append("\t").append((double) gerpScore[i] / gerpAlignCount[i]).append("\t");
                }

                if (snpGerpAlignCount[i] == 0) {
                    sb.append(Double.NaN).append("\t").append(Double.NaN);
                } else {
                    sb.append((double) snpGerpTree[i] / snpGerpAlignCount[i]).append("\t").append((double) snpGerpScore[i] / snpGerpAlignCount[i]);
                }

                double cdsL = (double) (snpCount[i] - noAncCount[i]) / snpCount[i] * cdsLength;
                sb.append("\t").append(noAncCount[i]).append("\t").append(b73SynCount[i]).append("\t").append((double) b73SynCount[i] / cdsL).append("\t");
                sb.append(b73NonCount[i]).append("\t").append((double) b73NonCount[i] / cdsL).append("\t");

                sb.append(b73DelCount[i]).append("\t").append((double) b73DelCount[i] / cdsL).append("\t");
                sb.append(b73DelHGCount[i]).append("\t").append((double) b73DelHGCount[i] / cdsL);

                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void scriptAddAncAllele() {
//        for (int i = 1; i < 43; i++) {
//            String CHR = PStringUtils.getNDigitNumber(3, i);
//            System.out.println("java -Xms50g -Xmx100g -jar 017_mkAnnoDB.addAncAllele.single.jar /data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/chr" 
//                    + CHR + ".lineage.maf0.005.bi.AnnoDB.txt.gz "
//                    + "/data4/home/aoyue/vmap2/analysis/ancestralAllele/Chr" + CHR 
//                    + ".ancestralAllele.txt " 
//                    + "/data4/home/aoyue/vmap2/analysis/015_annoDB/002_addAncestralAllele/chr"
//                    + CHR + ".lineage.maf0.005.bi.AnnoDB.addAncAllele.txt.gz"
//                    );
//        }

        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(l, i);
            if (index > -1) { //如果大于-1，则在集合中搜索到，说明是属于D基因组的；
                System.out.println("java -Xms50g -Xmx100g -jar 017_mkAnnoDB.addAncAllele.single.jar /data4/home/aoyue/vmap2/analysis/015_annoDB/003_addSIFT/chr"
                        + chr + ".lineage.maf0.005.bi.AnnoDB.addSIFT.txt.gz "
                        + "/data4/home/aoyue/vmap2/daxing/ancestralAllele/chr" + chr
                        + ".wheat.ancestralAllele.txt "
                        + "/data4/home/aoyue/vmap2/analysis/015_annoDB/005_addAncestralAllele/chr"
                        + chr + ".lineage.maf0.005.bi.AnnoDB.addSIFT.addAnc.txt.gz > log_017/log_" + chr + "_mkAnnoDB.addAncAllele.single.txt"
                );
            }
        }
    }

    /**
     *
     * @param dbfileS
     * @param ancS
     * @param outfileS
     */
    public void addAncAllele_singlethread(String dbfileS, String ancS, String outfileS) {
        double daf = Double.NaN;

        //Step1:建立42个文件的snpPos的二维数组，第一维是染色体号，第二维是每条染色体含有的pos集合；
        //Step2:然后文件并行流读进去，将染色体号提取出来，提取的chr的数字形式-1就是snpPos[][]的第一维；
        //在文件流内建立一个集合posList，将pos信息加入posList中，文件读入完毕后，将该posList转换成数组，放入每一维染色体对应的pos信息。
        //在文件流内建立一个HashMap集合，将pos对应的ancestral allele建立联系，
        //Step3:读入数据库文件，进行一行一行搜索，如果posList中有该位点，说明存在ancestral alle,即在文件最后一列添加，如果没有搜到，则输出NA
        File f = new File(ancS);
        int chr = Integer.parseInt(f.getName().substring(3, 6));
        TIntArrayList snpPosList = new TIntArrayList();
        HashMap<Integer, String> hm = new HashMap<>();

        try { //Chr001.ancestralAllele.txt  chr001.lineage.maf0.005.bi.AnnoDB.txt.gz
            BufferedReader br = null;
            if (f.getName().endsWith(".txt")) {
                br = IOUtils.getTextReader(ancS);
            } else if (f.getName().endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(ancS);
            }
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String anc = l.get(3);
                snpPosList.add(pos);
                hm.put(pos, anc);
            }
            br.close();
            System.out.println(f.getName() + "\tis completed on posList DB ");
        } catch (Exception e) {
            e.printStackTrace();
        }

        int[] snpPos = snpPosList.toArray(new int[snpPosList.size()]);
        Arrays.sort(snpPos);

        try { //chr001.lineage.maf0.005.bi.AnnoDB.txt.gz
            String chrS = PStringUtils.getNDigitNumber(3, chr);
            BufferedReader br = null;
            if (dbfileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(dbfileS);
            } else if (dbfileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(dbfileS);
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //压缩格式的输出
            String temp = br.readLine();
            bw.write(temp + "\tAncestralAllele\tDaf");
            bw.newLine();
            int cntAnc = 0;
            int cntAncNotMajororMinor = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                String major = l.get(4);
                String minor = l.get(5);
                double maf = Double.parseDouble(l.get(6));
                int index = Arrays.binarySearch(snpPos, pos);
                StringBuilder sb = new StringBuilder();
                if (index > -1) { //表明含有anc
                    String ancAllele = hm.get(pos);
                    //如果ancestral allele存在,且等于major，则da等于的minor, daf 就等于maf
                    //如果ancestral allele存在,且等于minor，则da等于的major, daf 就等于 1-daf1
                    if (ancAllele.equals(minor)) {
                        cntAnc++;
                        daf = 1 - maf;
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (ancAllele.equals(major)) {
                        cntAnc++;
                        daf = maf;
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (!ancAllele.equals(minor) && !ancAllele.equals(major)) {
                        sb.append(temp).append("\t").append(ancAllele).append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                        cntAncNotMajororMinor++;
                    }

                } else { //表明不含anc
                    sb.append(temp).append("\t").append("NA").append("\t").append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
            System.out.println("chr" + PStringUtils.getNDigitNumber(3, chr) + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor.");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void addAncestralAllele(String infileDirS, String ancestralDirS, String outfileDirS) {
        //String infileDirS = "";
        //String outfileDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/002_addAncestralAllele";

        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }

        //Step1:建立42个文件的snpPos的二维数组，第一维是染色体号，第二维是每条染色体含有的pos集合；
        //Step2:然后文件并行流读进去，将染色体号提取出来，提取的chr的数字形式-1就是snpPos[][]的第一维；
        //在文件流内建立一个集合posList，将pos信息加入posList中，文件读入完毕后，将该posList转换成数组，放入每一维染色体对应的pos信息。
        //在文件流内建立一个HashMap集合，将pos对应的ancestral allele建立联系，
        //Step3:读入数据库文件，进行一行一行搜索，如果posList中有该位点，说明存在ancestral alle,即在文件最后一列添加，如果没有搜到，则输出NA
        File[] fs = new File(ancestralDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(ancestralDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            int chr = Integer.parseInt(f.getName().substring(3, 6));
            TIntArrayList snpPosList = new TIntArrayList();
            HashMap<Integer, String> hm = new HashMap<>();

            try { //Chr001.ancestralAllele.txt  chr001.Alineage.maf0.005.bi.AnnoDB.txt.gz
                String infileS = f.getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = br.readLine(); //read header
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String anc = l.get(3);
                    snpPosList.add(pos);
                    hm.put(pos, anc);
                }
                br.close();
                System.out.println(f.getName() + "\tis completed on posList DB ");
            } catch (Exception e) {
                e.printStackTrace();
            }

            int[] snpPos = snpPosList.toArray(new int[snpPosList.size()]);
            Arrays.sort(snpPos);

            try { //chr001.Alineage.maf0.005.bi.AnnoDB.txt.gz
                String chrS = PStringUtils.getNDigitNumber(3, chr);
                String infileS = new File(infileDirS, "chr" + String.valueOf(chrS) + "." + hml.get(chr) + "lineage.maf0.005.bi.AnnoDB.txt.gz").getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, "chr" + String.valueOf(chrS) + "." + hml.get(chr) + "lineage.maf0.005.bi.AnnoDB.addAncAllele.txt.gz").getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, "chr" + String.valueOf(chrS) + "." + hml.get(chr) + "lineage.maf0.005.bi.AnnoDB.addAncAllele.txt.gz").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //压缩格式的输出
                String temp = br.readLine();
                bw.write(temp + "\tAncestralAllele\tDaf");
                bw.newLine();
                int cntAnc = 0;
                int cntAncNotMajororMinor = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String major = l.get(4);
                    String minor = l.get(5);
                    double maf = Double.parseDouble(l.get(6));
                    int index = Arrays.binarySearch(snpPos, pos);
                    if (index > -1) { //表明含有anc
                        String ancAllele = hm.get(pos);
                        double daf = Double.NaN;
                        //如果ancestral allele存在,且等于major，则da等于的minor, daf 就等于maf
                        //如果ancestral allele存在,且等于minor，则da等于的major, daf 就等于 1-daf1
                        if (ancAllele.equals(minor)) {
                            cntAnc++;
                            daf = 1 - maf;
                            StringBuilder sb = new StringBuilder();
                            sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                        if (ancAllele.equals(major)) {
                            cntAnc++;
                            daf = maf;
                            StringBuilder sb = new StringBuilder();
                            sb.append(temp).append("\t").append(ancAllele).append("\t").append(String.format("%.4f", daf));
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                        if (!ancAllele.equals(minor) && !ancAllele.equals(major)) {
                            StringBuilder sb = new StringBuilder();
                            sb.append(temp).append("\t").append(ancAllele).append("\t").append("NA");
                            bw.write(sb.toString());
                            bw.newLine();
                            //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                            // cntAncNotMajororMinor

                        }

                    } else if (index < 0) { //表明不含anc
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp).append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
                System.out.println("chr" + PStringUtils.getNDigitNumber(3, chr) + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor.");
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     *
     *
     * @param infileDirS
     */
    public void mkSNPsummary(String infileDirS, String outfileDirS) {

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
        /**
         * *************************************************************
         */
        System.out.println("FileName\tbiallelicNum\tBiallelicMafmore0.005Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".AnnoDB.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".AnnoDB.txt.gz")).getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tHetProportion\tMissProportion");
                bw.newLine();
                String temp = null;
                String te[] = null;
                int biallelicNum = 0;
                int biallelicMafmoreNum = 0;

                String major = null;
                String minor = null;
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
                        //bw.write(temp);
                        //bw.newLine();
                    } else {
                        te = temp.split("\t");
                        String chr = PStringUtils.fastSplit(temp).get(0);
                        String pos = PStringUtils.fastSplit(temp).get(1);
                        String ref = PStringUtils.fastSplit(temp).get(3);
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
                            hetRate = hetNum / genoNum;
                            missRate = missNum / (missNum + genoNum);
                            refAF = refAlleleGametes / (refAlleleGametes + altAlleleGametes);
                            altAF = altAlleleGametes / (refAlleleGametes + altAlleleGametes);;
                            if (refAF >= altAF) {
                                major = ref;
                                minor = alt;
                                maf = altAF;
                            } else {
                                maf = refAF;
                                major = alt;
                                minor = ref;
                            }

                            if (maf <= 0.005) {
                                continue;
                            }
                            biallelicMafmoreNum++;
                            StringBuilder sb = new StringBuilder();
                            //bw.write("Chr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tHetProportion\tMissProportion");
                            sb.append(chr).append("\t").append(pos).append("\t").append(ref).append("\t").append(alt).append("\t").
                                    append(major).append("\t").append(minor).append("\t").append(String.format("%.4f", maf)).append("\t").
                                    append(String.format("%.4f", hetRate)).append("\t").append(String.format("%.4f", missRate));
                            bw.write(sb.toString());
                            bw.newLine();

                        }
                    } //else的终止
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

}
