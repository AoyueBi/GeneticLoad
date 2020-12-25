package AoUtils;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.table.RowTable;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author AoyueBi
 * @data 2020-10-21 16:47
 */
public class AoGene {

    public AoGene(){

    }

    /**
     * 根据染色体的位置来返回对应的基因
     * @param chrIDArray
     * @param posArray
     */
    public static String[] getTranscriptName(int[] chrIDArray, int[] posArray) {
        String[] out = new String[chrIDArray.length];
        List<String> outtransList = new ArrayList<>();
        String geneHCFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        RowTable<String> t = new RowTable<>(geneHCFileS);

        //// build list array and initialize it
        int chrNum = 42;
        List<String>[] geneList = new ArrayList[chrNum];
        List<String>[] tranList = new ArrayList[chrNum];
        TIntList[] startLists = new TIntList[chrNum]; //所有的起始位点建立一个集合
        TIntList[] endLists = new TIntList[chrNum]; //所有的终止位点建立一个集合

        for (int i = 0; i < chrNum; i++) { //对list数组进行初始化
            geneList[i] = new ArrayList<>();
            tranList[i] = new ArrayList<>();
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
        }

        //// add value to the list and finish the database building
        int chrIndex = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int currentChr = Integer.parseInt(t.getCell(i, 2));
            chrIndex = currentChr - 1;
            geneList[chrIndex].add(t.getCell(i, 0));
            tranList[chrIndex].add(t.getCell(i, 1));
            startLists[chrIndex].add(Integer.parseInt(t.getCell(i, 3))); //本条染色体内的所有基因的起始位置组合
            endLists[chrIndex].add(Integer.parseInt(t.getCell(i, 4)));
        }

        for (int i = 0; i < startLists.length; i++) {
            startLists[i].sort();
            endLists[i].sort();
        }

        System.out.println("Finished building the gene repositry list");


        //// start to analysis the gene name of a query pos
        for (int i = 0; i < chrIDArray.length; i++) {
            int index = chrIDArray[i] - 1; //染色体号的索引
            int pos = posArray[i];
            int posIndex = -1;
            posIndex = startLists[index].binarySearch(pos);
            if (posIndex < 0) {
                posIndex = -posIndex - 2; //确保该位点在起始位点的右边
            }
            if (posIndex < 0) {
                outtransList.add("NA");
                continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
            }
            if (pos >= endLists[index].get(posIndex)) {
                outtransList.add("NA");
                continue; //确保在末端位点的前面，若不在，也舍去
            }
            outtransList.add(tranList[index].get(posIndex));
        }

        out = outtransList.toArray(new String[outtransList.size()]);

        return out;
    }

}

