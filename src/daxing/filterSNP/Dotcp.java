package daxing.filterSNP;

import pgl.infra.pos.ChrPos;

public class Dotcp extends ChrPos {

    private double depth;
    private double sd;

    public Dotcp(short chr, int pos, double depth, double sd){
        super(chr, pos);
        this.depth=depth;
        this.sd = sd;
    }

    public double getDepth(){
        return depth;
    }

    public double getSd() {
        return sd;
    }
}
