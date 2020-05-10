package WheatGeneticLoad;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import com.koloboke.collect.map.hash.HashIntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;

import java.util.Arrays;

class Grid {
    double xStart;
    double xEnd;
    double yStart;
    double yEnd;
    int nBin;

    double[] xBound;
    double[] yBound;
    int[] coordinates;
    int[] cellCounts;
    HashIntIntMap coordinateOrderMap = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();

    public Grid (double xStart, double xEnd, double yStart, double yEnd, int nBin) {
        this.xStart = xStart;
        this.xEnd = xEnd;
        this.yStart = yStart;
        this.yEnd = yEnd;
        this.nBin = nBin;
        xBound = this.initializeBin(xStart, xEnd);
        yBound = this.initializeBin(yStart, yEnd);
        this.initilializeCell();
    }

    private void initilializeCell () {
        this.coordinates = new int[nBin*nBin];
        this.cellCounts = new int[nBin*nBin];
        int cnt = 0;
        for (int i = 0; i < nBin; i++) {
            for (int j = 0; j < nBin; j++) {
                this.coordinates[cnt] = this.getCoordinate((short)i, (short)j);
                cnt++;
            }
        }
    }

    private int getCoordinate (short x, short y) {
        int value = x << 16;
        value = value+y;
        return value;
    }

    private short[] decodeCoordinate (int coordinate) {
        int x = coordinate >> 16;
        int v = ~0;
        v = v >> 16;
        int y = v & coordinate;
        short[] c = new short[2];
        c[0] = (short)x;
        c[1] = (short)y;
        return c;
    }

    private double[] initializeBin (double start, double end) {
        double interval = (end-start)/nBin;
        double[] bound = new double[nBin+1];
        for (int i = 0; i < bound.length; i++) {
            bound[i] = start+i*interval;
        }
        return bound;
    }

    public int getTotalSiteCount () {
        int cnt = 0;
        for (int i = 0; i < this.cellCounts.length; i++) {
            cnt+=cellCounts[i];
        }
        return cnt;
    }

    public int getOrderIndexOfProportionOfSite (double proportion) {
        int totalCount = this.getTotalSiteCount();
        double current = 0;
        for (int i = cellCounts.length-1; i > -1; i--) {
            current+=(double)cellCounts[i]/totalCount;
            if (current>proportion) return i;
        }
        return Integer.MIN_VALUE;
    }

    public boolean isHighDensity (double x, double y, double indexThresh) {
        int xIndex = Arrays.binarySearch(this.xBound, x);
        if (xIndex < 0) {
            xIndex = -xIndex-2;
        }
        if (xIndex < 0) return false;
        if (xIndex > nBin - 1) return false;
        int yIndex = Arrays.binarySearch(this.yBound, y);
        if (yIndex < 0) {
            yIndex = -yIndex-2;
        }
        if (yIndex < 0) return false;
        if (yIndex > nBin - 1) return false;
        int coor = this.getCoordinate((short)xIndex, (short)yIndex);
        int index = this.coordinateOrderMap.get(coor);
        if (index > indexThresh) return true;
        return false;
    }

    public void addXY (double x, double y) {
        int xIndex = Arrays.binarySearch(this.xBound, x);
        if (xIndex < 0) {
            xIndex = -xIndex-2;
        }
        if (xIndex < 0) return;
        if (xIndex > nBin - 1) return;
        int yIndex = Arrays.binarySearch(this.yBound, y);
        if (yIndex < 0) {
            yIndex = -yIndex-2;
        }
        if (yIndex < 0) return;
        if (yIndex > nBin - 1) return;
        int index = xIndex*nBin+yIndex;
        this.cellCounts[index]++;
    }

    private void sortByCellcounts () {
        GenericSorting.quickSort(0, cellCounts.length, compByCount, swapper);
    }

    public void buildHashMap () {
        this.sortByCellcounts();
        for (int i = 0; i < nBin; i++) {
            for (int j = 0; j < nBin; j++) {
                int coor = this.getCoordinate((short)i, (short)j);
                for (int k = 0; k < coordinates.length; k++) {
                    if (coor != coordinates[k]) continue;
                    this.coordinateOrderMap.put(coor, k);
                    break;
                }

            }
        }
    }

    protected Swapper swapper = new Swapper() {
        @Override
        public void swap(int a, int b) {
            int temp = cellCounts[a];
            cellCounts[a] = cellCounts[b];
            cellCounts[b] = temp;
            temp = coordinates[a];
            coordinates[a] = coordinates[b];
            coordinates[b] = temp;
        }
    };

    protected IntComparator compByCount = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            if (cellCounts[a] < cellCounts[b]) return -1;
            else if (cellCounts[a] > cellCounts[b]) return 1;
            else return 0;
        }
    };
}


