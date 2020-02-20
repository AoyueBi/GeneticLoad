package Annotation;

import analysis.wheatVMap2.VMapDBUtils;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.format.table.ColumnTable;
import pgl.format.table.RowTable;
import pgl.utils.Dyad;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;
import pgl.utils.wheat.RefV1Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;

public class AnnotationCrossover {
    public AnnotationCrossover(){
        this.convertCoordinate();

    }

    /**
     *
     *
     */
    public void addRecombination () {
        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/test/001";
        String outDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/test/002_addRecombination";
        String recombinationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_recombination_rate_chrID.txt";
        ColumnTable<String> t = new ColumnTable<>(recombinationFileS);
        int chrNum = Integer.parseInt(t.getCell(t.getRowNumber()-1, 0));
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

    /**
     * 通过 Chr1A- Chr7D pos 信息转换为 chrID 格式
     *
     */
    public void convertCoordinate () {
        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_recombination_rate.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_crossover_rate_chrID.txt";
        RowTable<String> t = new RowTable<>(infileS);
        String header = "ChrID\tPosStart\tPosEnd\tCrossover_frequency";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            String chromosome = null;
            int start = -1;
            int end = -1;
            int chr1 = -1;
            int chr2 = -1;
            int pos1 = -1;
            int pos2 = -1;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                chromosome = t.getCell(i,0).replaceFirst("chr", "");
                start = Integer.parseInt(t.getCell(i,1));
                end = Integer.parseInt(t.getCell(i,2));
                chr1 = RefV1Utils.getChrID(chromosome, start);
                chr2 = RefV1Utils.getChrID(chromosome, end);
                if (chr1 == chr2) { //如果所在区间的起始和末尾都属于一条染色体
                    pos1 = RefV1Utils.getPosOnChrID(chromosome, start);
                    pos2 = RefV1Utils.getPosOnChrID(chromosome, end);
                    sb.append(chr1).append("\t").append(pos1).append("\t").append(pos2).append("\t").append(t.getCell(i, 4));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                else {
                    pos1 = RefV1Utils.getPosOnChrID(chromosome, start);
                    sb.append(chr1).append("\t").append(pos1).append("\t").append(RefV1Utils.getChrIDLength(chr1)+1).append("\t").append(t.getCell(i, 4));
                    bw.write(sb.toString());
                    bw.newLine();
                    sb.setLength(0);
                    pos2 = RefV1Utils.getPosOnChrID(chromosome, end);
                    sb.append(chr2).append("\t").append(1).append("\t").append(pos2).append("\t").append(t.getCell(i, 4));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
