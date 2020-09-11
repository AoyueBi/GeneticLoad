package AoUtils.AoClass;

import WheatGeneticLoad.FilterVCF2;

/**
 * @author AoyueBi
 * @data 2020-09-10 00:41
 */
public class AoRecord implements Comparable<AoRecord> {
    int pos;
    public String r;
    public AoRecord (int pos, String r){
        this.pos = pos;
        this.r = r;
    }
    @Override
    public int compareTo(AoRecord o){return this.pos-o.pos;}

//
//    class Record implements Comparable<Record>{
//        int pos;
//        public String r;
//        public Record (int pos, String r){
//            this.pos = pos;
//            this.r = r;
//        }
//        @Override
//        public int compareTo(Record o){return this.pos-o.pos;}
//    }

}
