package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.CalVCF;
import AoUtils.WheatUtils;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import pgl.infra.window.SimpleWindow;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author AoyueBi
 * @data 2020-09-09 21:09
 */
public class FdVSdel {

    public FdVSdel(){

        //*********** Fd vs del ********************//
//        this.testFTTvsLR();
//        this.mergeFd();

//        this.getsubspeciesSNPAnnotation();
//        this.mergeExonSNPAnnotation();

        this.WindowDelvsSyn_fromExonAnnotation();

    }



    public void mergeExonSNPAnnotation(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/003_SNPannotation_LR_EU/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/003_SNPannotation_LR_EU/002/001_exonSNP_anno_landraceEU.txt.gz";
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    /**
     * 这里只研究 del 和 syn
     */
    public void WindowDelvsSyn_fromExonAnnotation(){
//        String infileS = "";
//        String outfileS = "";
        int windowSize = 2000000; //2 M
        int windowStep = 1000000; //1 M

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/003_SNPannotation_LR_EU/002/001_exonSNP_anno_landraceEU.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/004_genomeScan/001/001_delVSsynOnchr_" + windowSize + "Window_" + windowStep + "step.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/016_genomeScan_delvsSyn/001/001_delVSsynOnChr_" + windowSize + "Window" + windowStep + "step.txt";
        String outfile2S = outfileS.replaceFirst(".txt","_addEffectiveCDSLength.txt");
        AoFile.readheader(infileS);
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
//        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};

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
//*********************************** add recombination *******************************************
            this.addCDSLengthInWindow(outfileS, windowSize, windowStep,outfile2S); //*********************************** new method *******************************
            System.out.println("======== " + "add effective CDS length.");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
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
            for (int j = 0; j < gt.getRowNumber(); j++) { //
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
            } //for
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
     * 根据亚群的个体，结合merged VCF 文件，提取亚群的 snp annotation info
     *
     */
    public void getsubspeciesSNPAnnotation(){
//        String taxaList = "";
//        String outfileDirS = "";
        String vcfFileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";


        //需要修改
        //landrace_EU
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/012_GroupbyContinent/LR_EU.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/003_SNPannotation_LR_EU/001";

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
            String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_landraceEU_.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_cultivar_.txt.gz").getAbsolutePath();
            int colmnIndex = 2;
            File out = AoFile.filterTxtLines(infileS,2,posl,outfileS);
        });

    }

    public void mergeFd(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/002_Fd_FTTvsLR";
        String outfileS= "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/002/Fd_vmap2.1.geno.Landrace.Free_threshing_tetraploid.txt";
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    /**
     * 获取ftt -> 欧洲LR 渗入的结果
     */
    public void testFTTvsLR(){ //LR指的是欧洲的LR
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/001_fdBySubspeciesPop_2MWindow1MStep";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/040_fdVSdel/002_Fd_FTTvsLR";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fsList.size(); i++) {
            File f = fsList.get(i);
            String filename = f.getName();
            if (!filename.contains("Landrace.Free_threshing_tetraploid") && (!filename.contains("Landrace.Ae.tauschii")))continue;
            fList.add(f);
        }

        fList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName()).getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);bw.newLine();
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

}
