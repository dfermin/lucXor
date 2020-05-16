package lucxor;

import java.util.*;

/**
 * PSM List that is use by the algorithm to
 */
public class PSMList implements List<PSM>{

  // List of PSMs
  List<PSM> psmList = new ArrayList<>();

  // This hash map controls the FileName + Scan, Order of the PSM in the psmList
  Map<String, Integer> scanOrder = new HashMap<>();

  int numDecoys = 0;
  int numTargets = 0;

  public int getNumDecoys() {
    return numDecoys;
  }

  public int getNumTargets() {
    return numTargets;
  }

  @Override
  public int size() {
    return psmList.size();
  }

  @Override
  public boolean isEmpty() {
    return psmList.isEmpty();
  }

  @Override
  public boolean contains(Object o) {
    return psmList.contains(o);
  }

  @Override
  public Iterator<PSM> iterator() {
    return psmList.iterator();
  }

  @Override
  public Object[] toArray() {
    return psmList.toArray();
  }

  @Override
  public <T> T[] toArray(T[] a) {
    return psmList.toArray(a);
  }


  @Override
  public boolean add(PSM psm) {
    int index = psmList.size();
    psmList.add(psm);
    scanOrder.put(Utils.generateIndex(psm.getSrcFile(), psm.getScanNum()), index);
    if (psm.isDecoy())
      numDecoys++;
    else
      numTargets++;
    return true;
  }

  @Override
  public boolean remove(Object o) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public boolean addAll(Collection c) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public boolean addAll(int index, Collection c) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public void clear() {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public boolean retainAll(Collection c) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public boolean removeAll(Collection c) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public boolean containsAll(Collection c) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public PSM get(int index) {
    return psmList.get(index);
  }

  @Override
  public PSM set(int index, PSM element) {
    throw new UnsupportedOperationException("Non supported metho");
  }

  @Override
  public void add(int index, PSM element) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public PSM remove(int index) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public int indexOf(Object o) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public int lastIndexOf(Object o) {
    throw new UnsupportedOperationException("Non supported method");
  }

  @Override
  public ListIterator<PSM> listIterator() {
    return psmList.listIterator();
  }

  @Override
  public ListIterator<PSM> listIterator(int index) {
    return psmList.listIterator(index);
  }

  @Override
  public List<PSM> subList(int fromIndex, int toIndex) {
    return psmList.subList(fromIndex, toIndex);
  }

  /**
   * This function retrieve a PSM by the scanNumber in the File
   * @param baseFN FileName
   * @param scanNum Scan Number
   * @return PSM
   */
  public PSM getByScanOrder(String baseFN, int scanNum) {
    String hash = Utils.generateIndex(baseFN, scanNum);
    int index = scanOrder.get(hash);
    return psmList.get(index);
  }
}