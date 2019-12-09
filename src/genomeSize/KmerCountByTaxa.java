/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package genomeSize;

import com.koloboke.collect.set.hash.HashLongSet;
import com.koloboke.collect.set.hash.HashLongSets;
import format.dna.BaseEncoder;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class KmerCountByTaxa {

    String[] taxa = null;
    int[][] count = null;
    List<int[]> ls = null;
    long[] kmers = null;
    int kmerLength = 0;

    public KmerCountByTaxa(KmerCount[] kcs) {
        this.initialize(kcs);
    }

    public void initialize(KmerCount[] kcs) {
        taxa = new String[kcs.length];
        for (int i = 0; i < taxa.length; i++) {
            taxa[i] = kcs[i].getTaxon();
        }
        HashLongSet kmerSet = HashLongSets.newMutableSet();
        for (int i = 0; i < kcs.length; i++) {
            kmerSet.addAll(HashLongSets.newImmutableSet(kcs[i].kmers));
        }
        kmers = kmerSet.toLongArray();
        
        int index = 0;
        for (int i = 0; i < kmers.length; i++) {
            for (int j = 0; j < kcs.length; j++) {
                index = kcs[j].getKmerIndex(kmers[i]);
                if (index < 0) {
                    continue;
                }
                int a = kcs[j].getKmerCount(index);
                ls.get(i)[j] = a;
            }
        }
        kmerLength = kcs[0].getKmerLength();
    }

//    public void initialize (KmerCount[] kcs) {
//        taxa = new String[kcs.length];
//        for (int i = 0; i < taxa.length; i++)         {
//            taxa[i] = kcs[i].getTaxon();
//        }
//        HashLongSet kmerSet = HashLongSets.newMutableSet();
//        for (int i = 0; i < kcs.length; i++) {
//            kmerSet.addAll(HashLongSets.newImmutableSet(kcs[i].kmers));
//        }
//        kmers = kmerSet.toLongArray();
//        count = new int[kcs.length][kmers.length];
//        int index = 0;
//        for (int i = 0; i < kmers.length; i++) {
//            for (int j = 0; j < kcs.length; j++) {
//                index = kcs[j].getKmerIndex(kmers[i]);
//                if (index < 0) continue;
//                count[i][j] = kcs[j].getKmerCount(index);
//            }
//        }
//        kmerLength = kcs[0].getKmerLength();
//    }
    
    
    public void writeTextFile(String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Kmers");
            for (int i = 0; i < taxa.length; i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                sb.setLength(0);
                sb.append(BaseEncoder.getSequenceFromLong(kmers[i]).substring(0, this.kmerLength));
                for (int j = 0; j < count.length; j++) {
                    sb.append("\t").append(count[i][j]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
//    public void writeTextFile(String outfileS) {
//        try {
//            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//            StringBuilder sb = new StringBuilder("Kmers");
//            for (int i = 0; i < taxa.length; i++) {
//                sb.append("\t").append(taxa[i]);
//            }
//            bw.write(sb.toString());
//            bw.newLine();
//            for (int i = 0; i < count.length; i++) {
//                sb.setLength(0);
//                sb.append(BaseEncoder.getSequenceFromLong(kmers[i]).substring(0, this.kmerLength));
//                for (int j = 0; j < count.length; j++) {
//                    sb.append("\t").append(count[i][j]);
//                }
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            bw.flush();
//            bw.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
}
