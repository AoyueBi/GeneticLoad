package WheatGeneticLoad;

import AoUtils.AoFile;
import analysis.wheatVMap2.DBGene;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.sqlite.core.DB;
import pgl.format.genomeAnnotation.GeneFeature;
import pgl.format.range.Range;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class DBgene {

    public DBgene(){
        this.getTranscriptSum();

    }

    /**
     * 获取转录本的总结
     */
    public void getTranscriptSum(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/001/transcriptSummary.txt";

        double gerpCut = 0;
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
//        for (int i = 0; i < gf.getGeneNumber(); i++) {
//            int chrIndex = gf.getGeneChromosome(i)-1;
//            if (! (chrIndex== 0))continue;
//            int longTransIndex = gf.getLongestTranscriptIndex(i); //根据基因获取最长转录本索引
//            String geneName = gf.getTranscriptName(i, longTransIndex); //根据索引获取最长转录本的名字
//            genesList.add(geneName); //将最长转录本的名字加入geneList
//            List<Range> cdsList = gf.getCDSList(i, longTransIndex); //根据索引获取最长转录本的CDSList
//            int cnt = 0;
//
//            /*对于每一个基因的编码序列，还有很多个cds片段，即cdsList；我们对cdsList进行for循环，得到每个cds的起始和终止位置，从而计算出总长*/
//            for (int j = 0; j < cdsList.size(); j++) {
//                int rStart = cdsList.get(j).start;
//                int rEnd = cdsList.get(j).end;
//                for (int k = rStart; k < rEnd; k++) {
//                    /*posGeneMap是一个HashMap数组，一条染色体对应一个String类型的ArrayList；
//                    故得到该位点的所属基因名字列表，如果该位点不含基因名，就将 genename赋值给该位点，完善posGeneMap
//                    否则，如果该位点含有其他基因的基因名字，依旧把genename赋值给该位点*/
//                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(k); //建立map的关系，那个位点对应哪个list HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
//                    if (geneNameList == null) {
//                        geneNameList = new ArrayList();
//                        geneNameList.add(geneName);
//                        posGeneMap[chrIndex].put(k, geneNameList);
//                    }
//                    else {
//                        geneNameList.add(geneName);
//                        posGeneMap[chrIndex].put(k, geneNameList); /*最终将posGeneMap绘图完成*/
//                    }
//                    cnt++; /*每一个CDS位点相加，最终得到这个cds的长度。*/
//                }
//
//                // 最终cnt是一个基因的所有cdslist中，每个cds的每个位点包含的基因数目的总和
//            } //该循环是一个基因的所有cds循环
//            geneCDSLengthMap.put(geneName, cnt); //
//        }

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

            for (int i = 0; i < gf.getGeneNumber(); i++) {
                int chrindex = gf.getGeneChromosome(i)-1;
                if (! (chrindex== chrIndex))continue;
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
}
