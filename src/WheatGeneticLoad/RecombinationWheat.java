/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import gnu.trove.list.array.TIntArrayList;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class RecombinationWheat {

    public RecombinationWheat() {

    }

    /**
     * 弃用
     *
     */
    private void mkRecombinationPointTable() {
        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_mapping_data.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_mkRecomninationPoint/recombinationPoint.txt";
        int chrNum = 21;
        TIntArrayList[] posList = new TIntArrayList[chrNum];

        for (int i = 0; i < posList.length; i++) { //每一个新的list都要new一下！！！因为上文只说明建立一个数组posList,大小为10.没有针对一个posList进行初始化。
            posList[i] = new TIntArrayList();
        }
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int chrIndex = Integer.valueOf(l.get(1));
                int pos = (Integer.valueOf(l.get(2)));
                posList[chrIndex].add(pos);
            }
            br.close();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos");
            bw.newLine();
            for (int i = 0; i < posList.length; i++) {
                int[] pos = posList[i].toArray();
                Arrays.sort(pos);
                for (int j = 0; j < pos.length; j++) {
                    bw.write(String.valueOf(i + 1) + "\t" + String.valueOf(pos[j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
