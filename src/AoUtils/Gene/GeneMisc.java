package AoUtils.Gene;

import AoUtils.AoFile;
import AoUtils.GeneFeatureAo;
import daxing.common.wheat.PGF;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GeneMisc {

    public GeneMisc(){
        this.getCDSRange();

    }



    /**
     * 按照bed 格式，输出 gene 的所有 cds 的 range 的 bed 格式
     * bed的每一行记录是左闭右开区间，并且以0作为染色体的起始；而gtf的每一行记录是包括区间的两端的，并且以1作为染色体的起始。
     */
    public void getCDSRange(){
        GeneFeatureAo gf = new GeneFeatureAo("/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf");
        gf.sortGeneByStartPosition();
        int cntGeneLen = 0;
        int cntCDSLen = 0;
        int geneNum = 0;
        //求所有基因的长度之和

        List<String> nonoverlapGeneList = AoFile.getStringList("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/006_nonoverlapGene/nonoverlapGene.txt",0);
        Collections.sort(nonoverlapGeneList);
        try{
            BufferedWriter bw = AoFile.writeFile("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/007_bed_onWholeGenome/GFF3_CDS_nonoverlap.bed");
            bw.write("Chr\tPos_start\tPos_end\tGeneName");
            bw.newLine();

            BufferedWriter bw2 = AoFile.writeFile("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/007_bed_onWholeGenome/GFF3_gene_nonoverlap.bed");
            bw2.write("Chr\tPos_start\tPos_end\tGeneName");
            bw2.newLine();

            for (int i = 0; i < gf.getGeneNumber(); i++) {
                int chr = gf.getChromosomeOfGene(i);
                if (chr==0)continue; ///过滤0号染色体的基因
                int j = gf.getLongestTranscriptIndex(i);
                String genename = gf.getGeneName(i);
                String trans = gf.getTranscriptName(i,j);
                int ifoverlap = Collections.binarySearch(nonoverlapGeneList,trans);
                if (ifoverlap < 0) continue; ///过滤重叠的基因

                int geneLen = gf.getGeneLength(i);
                int cdsLen = gf.getCDSLen(i,j);
                cntGeneLen += geneLen;
                cntCDSLen += cdsLen;
                geneNum++;

                int a1 = gf.getGeneStart(i);
                int a2 = gf.getGeneEnd(i); //注意 pgf 类 中的range 是左开右闭，故转换成 bed 格式，不需要进行 +1
                bw2.write(chr + "\t" + a1 + "\t" + a2 + "\t" + genename);
                bw2.newLine();

                StringBuilder sb = new StringBuilder();
                List<Range> cdsList = gf.getCDSList(i,j);
                for (int k = 0; k < cdsList.size(); k++) {
                    sb.setLength(0);
                    short chrcds = cdsList.get(k).chr;
                    int pos1 = cdsList.get(k).start;
                    int pos2 = cdsList.get(k).end; //注意 pgf 类 中的range 是左开右闭，故转换成 bed 格式，不需要进行 +1
                    sb.append(chrcds).append("\t").append(pos1).append("\t").append(pos2).append("\t").append(genename);
                    bw.write(sb.toString());
                    bw.newLine();
                }

            }
            bw2.flush();bw2.close();
            bw.flush();bw.close();
            System.out.println("Total gene length is " + cntGeneLen );
            System.out.println("Total cds length is " + cntCDSLen);
            System.out.println("Total gene num is " + geneNum);
            gf.getChromosomeList();

        }catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * 想要知道小麦基因组中整个基因占多少bp,占全基因组的比例
     */
    public void getGeneandCDSLength(){
        GeneFeatureAo gf = new GeneFeatureAo("/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf");
        gf.sortGeneByStartPosition();
        int cntGeneLen = 0;
        int cntCDSLen = 0;
        int geneNum = 0;
        //求所有基因的长度之和

        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chr = gf.getChromosomeOfGene(i);
            if (chr==0)continue;
            int j = gf.getLongestTranscriptIndex(i);
            int geneLen = gf.getGeneLength(i);
            int cdsLen = gf.getCDSLen(i,j);
            cntGeneLen += geneLen;
            cntCDSLen += cdsLen;
            geneNum++;
        }
        System.out.println("Total gene length is " + cntGeneLen );
        System.out.println("Total cds length is " + cntCDSLen);
        System.out.println("Total gene num is " + geneNum);
        gf.getChromosomeList();
    }


    /**
     * 在wheat_v1.1_nonoverlap.txt文件中添加列信息：Chr,TransStart,TransEnd,TranStrand,CDSExonNumber,CDSLength六列信息
     * //列出所要建立数据库的基因的详细信息表格, 在基因表格中添加 cds 长度，基因起始位点，exon个数等。
     */
    public void geneInfo(){
        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        PGF pgf = new PGF(geneFeatureFileS);
        pgf.sortGeneByName();
        gf.sortGeneByName();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tChr\tTransStart\tTransEnd\tTranStrand\tExonNumber\tCDSLength");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(0);
                String trans = l.get(4);
                int geneIndex = gf.getGeneIndex(gene);
                int longIndex = gf.getLongestTranscriptIndex(geneIndex);
                String longTrans = gf.getTranscriptName(geneIndex,longIndex);
                if (!trans.equals(longTrans)) System.out.println(temp);
                int chr = gf.getChromosomeOfGene(geneIndex);
                if (chr==0)continue;
                int start = gf.getTranscriptStart(geneIndex,longIndex);
                int end = gf.getTranscriptEnd(geneIndex,longIndex);
                int strand = gf.getTranscriptStrand(geneIndex,longIndex);
                int exonNum = gf.getExonList(geneIndex,longIndex).size();
                int CDSlength = pgf.getCDSLen(geneIndex, longIndex);
                sb.append("\t").append(chr).append("\t").append(start).append("\t").append(end).append("\t").append(strand).append("\t").append(exonNum).append("\t").append(CDSlength);
                bw.write(temp + sb.toString());
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

}
