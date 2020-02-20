package Annotation;

import pgl.format.table.RowTable;
import pgl.utils.IOUtils;
import pgl.utils.wheat.RefV1Utils;

import java.io.BufferedWriter;

public class AnnotationCrossover {
    public AnnotationCrossover(){
        this.convertCoordinate();

    }

    /**
     *
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
                if (chr1 == chr2) { //如果所在区间
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
