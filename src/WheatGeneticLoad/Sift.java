/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import format.dna.FastaBit;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Sift {

    public Sift() {
        //this.addGeneBiotype();
        new FastaBit("/Users/Aoyue/Documents/chr036.fa.gz").getName(0);

    }
    
    
    /**
     * modify gtf file, add gene_biotype to the 9th colum
     */
    public void addGeneBiotype() {
        String InputFileS = "/Users/Aoyue/Documents/wheat_v1.1_Lulab.gtf";
        String OutputFileS = "/Users/Aoyue/Documents/wheat_v1.1_Lulab.addBiotype.gtf";
        try {
            BufferedReader br = IOUtils.getTextReader(InputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(OutputFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                //String plate = null;
                //plate = PStringUtils.fastSplit(temp).get(2);
                //if(plate.equalsIgnoreCase("exon") | plate.equalsIgnoreCase("CDS") | plate.equalsIgnoreCase("stop_codon") | plate.equalsIgnoreCase("start_codon")){
                List<String> l = null;
                l = PStringUtils.fastSplit(temp);
                StringBuilder sb = new StringBuilder();
                sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3)).append("\t").
                        append(l.get(4)).append("\t").append(l.get(5)).append("\t").
                        append(l.get(6)).append("\t").append(l.get(7)).append("\t").append(l.get(8)).append(" ")
                        .append("gene_biotype").append(" ").append("\"protein_coding\";");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
