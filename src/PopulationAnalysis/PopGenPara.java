/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PopulationAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class PopGenPara {

    public PopGenPara() {
        this.getTaxaSet();

    }

    //Advanced/improved cultivar
    //Breeding/Research Material
    //Cultivar
    //Landrace
    //Other
    //Traditional cultivar/Landrace
    //NA
    public void renameType() {
        

    }

    public void getTaxaSet() {
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/011_hapScanner/taxaRefBam/abd/002_taxaRefBam.ABDgenome.manual.addNAFU.txt";
            String outfileS = "/Users/Aoyue/Documents/BreadWheat.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            Set<String> s = new HashSet<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                l = PStringUtils.fastSplit(temp);
                s.add(l.get(0));
                System.out.println(l.get(0));
            }
            br.close();
            System.out.println(s.size() + " taxa");

            String[] sS = s.toArray(new String[s.size()]);
            Arrays.sort(sS);
            for (int i = 0; i < s.size(); i++) {
                bw.write(sS[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
