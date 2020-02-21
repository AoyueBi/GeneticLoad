package PopulationAnalysis;

import pgl.format.table.RowTable;
import pgl.utils.IOUtils;
import pgl.utils.wheat.RefV1Utils;

import java.io.BufferedWriter;

public class XPCLR {
    public XPCLR(){
        this.convertCoordinate();

    }

    /**
     * Today I  just read REF and not programming
     * Today I  just made annual summary and not programming
     */

    public void test(){

    }

    /**
     * 通过 Chr1A- Chr7D pos 信息转换为 chrID 格式
     * 输入文件是 genetic pos map
     *
     */
    public void convertCoordinate () {
        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_mapping_data.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_mapping_data_chrID.txt";
        RowTable<String> t = new RowTable<>(infileS);
        String header = "psId\tchromosome\tphysicalPosition\tgeneticPosition";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            String chromosome = null;
            String psID = null;
            String geneticPos = null;
            int pos = -1;
            int chrid = -1;
            int posid = -1;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                psID = t.getCell(i,0);
                chromosome = t.getCell(i,1).replaceFirst("chr", "");
                pos = Integer.parseInt(t.getCell(i,2));
                geneticPos = t.getCell(i,3);
                chrid = RefV1Utils.getChrID(chromosome, pos);
                posid = RefV1Utils.getPosOnChrID(chromosome, pos);
                sb.append(psID).append("\t").append(chrid).append("\t").append(posid).append("\t").append(geneticPos);
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
