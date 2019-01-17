#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<string>
using namespace std;

const int INI_WEIGHT=1000000;


/**
 * Element stored within the Binary Heap.
 */
typedef struct elt {

  /** user-defined information to be stored by id */
  int     id;            

  /** Key which represents the priority. */
  int     priority;
} ELEMENT, *ELEMENT_PTR;


class BinaryHeap {

 public:
  BinaryHeap (int);
  ~BinaryHeap ();

  bool isEmpty() { return (_n == 0); }
  int smallest();
  void insert (int, int); 
  void decreaseKey (int, int);

 private:
  int          _n;           // number of elements in binary heap
  ELEMENT_PTR  _elements;    // values in the heap
  int          *_pos;        // pos[i] is index into elements for ith value    

  long         _numComparisons;
  long         _numInsert;
  long         _numSwaps;
  long         _numSmallest;
  long         _numDecrease;

};

/** allocate heap of n elements */
BinaryHeap::BinaryHeap (int i) {
  _n = 0;  // initially none in the heap

  // simplify algorithm to consider position 1 as being the 'root'
  _elements = (ELEMENT *) calloc (i+1, sizeof (ELEMENT));
  _pos      = (int *) calloc (i+1, sizeof (int));
}

/** Destructor frees all internal storage. */
BinaryHeap::~BinaryHeap () {
  free (_elements);
  free (_pos);
  _n = -1;
}


/**
 * Return the vertex identifier of the smallest vertex in heap and 
 * readjust the heap.
 */
int BinaryHeap::smallest () {
  int id = _elements[1].id;
  int pIdx;

  //INC_SMALL;

  // heap will have one less entry, and we want to place leftover one
  // in proper location.
  ELEMENT last = _elements[_n];
  _n--;

  _elements[1] = last;

  pIdx = 1;
  int child = pIdx*2;
  while (child <= _n) {
    // select smaller of two children
    ELEMENT sm = _elements[child];
    if (child < _n) {
      //INC_COMP;
      if (sm.priority >  _elements[child+1].priority) {
	sm = _elements[++child];
      }
    }

    // are we in the right spot? Leave now
    //INC_COMP;
    if (last.priority <= sm.priority) { break; }

    // otherwise, swap and move up
    //INC_SWAP;
    _elements[pIdx] = sm;
    _pos[sm.id] = pIdx;

    pIdx = child;
    child = 2*pIdx;
  }

  // insert into spot vacated by moved element (or last one)
  _elements[pIdx] = last;
  _pos[last.id] = pIdx;
  return id;
}


void BinaryHeap::insert (int id, int priority) {
  int i;

  //INC_INSERT;

  // add to end of the heap. If 1 then the first element.
  i = ++_n;
  while (i > 1) {
    int       pIdx = i/2;
    ELEMENT   p    = _elements[pIdx];

    // are we in the right spot? Leave now
    //INC_COMP;
    if (priority > p.priority) { break; }

    // otherwise, swap and move up
    //INC_SWAP;
    _elements[i] = p;
    _pos[p.id] = i;
    i = pIdx;
  }

  // insert into spot vacated by moved element (or last one)
  _elements[i].id = id;
  _elements[i].priority = priority;
  _pos[id] = i;
}


void BinaryHeap::decreaseKey (int id, int newPriority) {
  int size = _n;

  //INC_DECREASE;

  // truncate heap (temporarily) and act like the binary heap up to
  // but not including this one is all that exists (cute, huh?) 
  _n = _pos[id] - 1;

  // now we insert and the binary heap is shuffled appropriately
  insert(id, newPriority);

  // since the newPriority must be lower, we can expand back and 
  // we still have a working binary heap
  _n = size;
}

class travelStep{
    public:
        string name;
        int enterTime,leaveTime;
};


class nodeInfo{
    public:
        string name;
        int happiness;
        int openTime, closeTime;
        int graphIdx;
};


class Graph{
private:
    int total_vertex;
    int total_edge;

    vector<bool> visited;
   
    void singleSourceMinPath(int source, vector<int> &Adist, vector<int> &APred);
    
    void gohome(int start, vector<int> &apath);

    void gohomeAns2(int start, int travelTime, vector<int> &apath, vector<bool> &staySation, vector<nodeInfo> nodeGroup, int startTime);
    
public:
    vector<vector<int>> adjmatrix;
    vector<int> pred, dist;
    Graph(int n):total_vertex(n){
        total_edge=0;
        adjmatrix.resize(total_vertex);
        visited.resize(total_vertex);
        pred.resize(total_vertex);
        dist.resize(total_vertex);
        for(int i=0;i< total_vertex; i++){
            adjmatrix[i].resize(total_vertex);
            pred[i]=-1;
            dist[i]=INI_WEIGHT;
            
        }
        for (int i=0; i<total_vertex;i++)
				for(int j=0;j<total_vertex; j++)
				adjmatrix[i][j]=INI_WEIGHT;


    }
    void travel(int start, int travelTime, vector<nodeInfo> nodeGroup, vector<int> &aPath);
    void travelAns2(int start, int travelTime, vector<nodeInfo> nodeGroup, vector<int> &aPath, vector<bool> &stayStation,int startTime);
    void InitPredDist();
    void AddEdge(int start, int end, int weight);
    void CalculateMinPath(int source ){
        singleSourceMinPath(source, dist, pred);
        //Print();
    }

    void InitVisited(){
        for(int i=0; i<visited.size(); i++)
            visited[i]=false;
    }

    void Print(){
        for(int i=0; i<dist.size(); i++){
            cout<<dist[i]<<" ";
        }
        cout<<endl;

        for(int i=0; i<pred.size(); i++){
            cout<<pred[i]<<" ";
        }
        cout<<endl;        

    }

    int CalculateHappyAns1(vector<nodeInfo> nodeGroup, vector<int> aPath);

    int CalculateHappyAns2(vector<nodeInfo> nodeGroup, vector<int> aPath, vector<bool> stayGroup, int startTime); 

    int CalculateTravelTimeAns1(vector<nodeInfo> nodeGroup, vector<travelStep> &stepGroup, vector<int> aPath, int startTime);
 
    int CalculateTravelTimeAns2(vector<nodeInfo> nodeGroup, vector<travelStep> &stepGroup, vector<int> aPath, int startTime, vector<bool> stayGroup);
};
int Graph::CalculateTravelTimeAns2(vector<nodeInfo> nodeGroup, vector<travelStep> &stepGroup, vector<int> aPath, int startTime, vector<bool> stayGroup){
    int costTime=0;
    int tmpNode,next,waitTime=0;
    travelStep tmpStep;
    for (int i=0; i<aPath.size(); i++){
        tmpNode=aPath[i];
        waitTime=0;

        if (stayGroup[i]==true){
            waitTime=nodeGroup[tmpNode].openTime-(costTime+startTime);
        }
        tmpStep.name=nodeGroup[tmpNode].name;
        tmpStep.enterTime= costTime+startTime;
        tmpStep.leaveTime= costTime+startTime+waitTime;
        stepGroup.push_back(tmpStep);

        costTime=costTime+waitTime;

        if ((i+1)<aPath.size()){
            next=aPath[i+1];
            costTime=costTime+adjmatrix[tmpNode][next];
        }
    }

    return costTime;
}


int Graph::CalculateTravelTimeAns1(vector<nodeInfo> nodeGroup, vector<travelStep> &stepGroup, vector<int> aPath, int startTime){
    int costTime=0;
    int tmpNode,next;
    travelStep tmpStep;
    for (int i=0; i<aPath.size(); i++){
        tmpNode=aPath[i];
        tmpStep.name=nodeGroup[tmpNode].name;
        tmpStep.enterTime= costTime+startTime;
        tmpStep.leaveTime= costTime+startTime;
        stepGroup.push_back(tmpStep);

        if ((i+1)<aPath.size()){
            next=aPath[i+1];
            costTime=costTime+adjmatrix[tmpNode][next];
        }
    }

    return costTime;

}

int Graph::CalculateHappyAns1(vector<nodeInfo> nodeGroup, vector<int> aPath){
    InitVisited();
    int totalHappy=0;
    int tmpNode=0;
    for (int i=0; i<aPath.size(); i++){
        tmpNode=aPath[i];
        if (!visited[tmpNode]){
            visited[tmpNode]=true;
            totalHappy=totalHappy+nodeGroup[tmpNode].happiness;
        }
    }
    return totalHappy;
}

int Graph::CalculateHappyAns2(vector<nodeInfo> nodeGroup, vector<int> aPath, vector<bool> stayGroup,int startTime){
    InitVisited();
    int totalHappy=0;
    int tmpNode=0;
    int next=0;
    for (int i=0; i<aPath.size(); i++){
        tmpNode=aPath[i];
        if (stayGroup[i]==true){
            startTime=nodeGroup[tmpNode].openTime;
        }
        if (!visited[tmpNode] && startTime>=nodeGroup[tmpNode].openTime && startTime<=nodeGroup[tmpNode].closeTime){
            visited[tmpNode]=true;
            totalHappy=totalHappy+nodeGroup[tmpNode].happiness;
        }
        if ((i+1)<aPath.size()){
            next=aPath[i+1];
            startTime=startTime+adjmatrix[tmpNode][next];
        }
    }
    return totalHappy;

}

void Graph::InitPredDist(){
    for(int i=0; i<total_vertex; i++){
        pred[i]=-1;
        dist[i]=INI_WEIGHT;
    }
}

void Graph::gohome(int start ,vector<int> & apath){
    int pre=pred[start];
    while(pre!=-1){      
        apath.emplace_back(pre);
        pre=pred[pre];
    }
    //cout<<"gohome"<<endl;
}

void Graph::gohomeAns2(int start, int travelTime, vector<int> &apath, vector<bool> &staySation, vector<nodeInfo> nodeGroup, int startTime){
    int tmpstart=start;
    int pre=pred[start];
    int tmpStartTime=startTime;
    int tmpReamin=travelTime;
    int waitTime=0;
    bool stay=false;
    while(pre!=-1){ 
        stay=false;     
        waitTime=0;
        apath.emplace_back(pre);

        if (tmpStartTime+adjmatrix[tmpstart][pre]<nodeGroup[pre].openTime ){
            waitTime=nodeGroup[pre].openTime-tmpStartTime-adjmatrix[tmpstart][pre];
            if (tmpReamin-adjmatrix[tmpstart][pre]-waitTime >= dist[pre] ){
                stay=true;
                tmpReamin=tmpReamin-adjmatrix[tmpstart][pre]-waitTime;
                tmpStartTime=nodeGroup[pre].openTime;
                                
            }
            else{
                stay=false;
                tmpReamin=tmpReamin-adjmatrix[tmpstart][pre];
                tmpStartTime=tmpStartTime+adjmatrix[tmpstart][pre];
            }

        }
        else{
            stay=false;
            tmpReamin=tmpReamin-adjmatrix[tmpstart][pre];
            tmpStartTime=tmpStartTime+adjmatrix[tmpstart][pre];
            
        }
            

        staySation.push_back(stay);
        tmpstart=pre;
        pre=pred[pre];
    }    
}

void Graph::AddEdge(int start, int end, int weight){
    total_edge++;
    adjmatrix[start][end]=weight;
	adjmatrix[end][start]=weight;  
}

void Graph::singleSourceMinPath(int source, vector<int> &Adist, vector<int> &APred){
    Adist[source]=0;
    BinaryHeap pq(total_vertex);

    for (int u=0; u< total_vertex; u++){
        pq.insert(u, Adist[u]);
    }
    int newLen=0;
    int adjNode=0;   
    int smallNode=0;
    while(!pq.isEmpty()){
        smallNode=pq.smallest();
        for (int adjNode=0;adjNode<total_vertex ;adjNode++){
            newLen=Adist[smallNode]+adjmatrix[smallNode][adjNode];
            if (newLen<Adist[adjNode]){
                pq.decreaseKey(adjNode, newLen);
                Adist[adjNode]=newLen;
                APred[adjNode]=smallNode;                
            }
        }        
    }    
}

void Graph::travel(int start, int travelTime, vector<nodeInfo> nodeGroup, vector<int> &aPath){
    vector<int> tmpSelect;
    vector<int>().swap(tmpSelect);
    visited[start]=true;
    for (int i=0; i<total_vertex; i++ ){
        if ((adjmatrix[start][i] !=INI_WEIGHT ) && (travelTime-adjmatrix[start][i] >= dist[i] ) && (visited[i]==false) ){
            if (tmpSelect.size()>0){
                if (nodeGroup[i].happiness> nodeGroup[tmpSelect[0]].happiness){                  
                    tmpSelect.push_back(i);
                    tmpSelect.erase(tmpSelect.begin());
                }
            }
            else{
                //cout<<"empty"<<endl;
                tmpSelect.push_back(i);
            }    
        }
    }

    if (tmpSelect.size()>0){
        int next=tmpSelect[0];
        int remain=travelTime-adjmatrix[start][next];
        visited[next]=true;
        aPath.push_back(next);
        travel(next, remain, nodeGroup, aPath);
    }
    else{
        gohome(start, aPath);
    }
}


void Graph::travelAns2(int start, int travelTime, vector<nodeInfo> nodeGroup, vector<int> &aPath, vector<bool> &stayStation, int startTime){
    vector<int> tmpSelect;
    vector<int>().swap(tmpSelect);
    visited[start]=true;
    int waitTime=0;
    bool findNext=false;
    for (int i=0; i<total_vertex; i++ ){
        if ((adjmatrix[start][i] !=INI_WEIGHT ) && (travelTime-adjmatrix[start][i] >= dist[i] ) && (visited[i]==false) ){
            if (tmpSelect.size()>0){
                if (startTime+adjmatrix[start][i]>= nodeGroup[i].openTime && 
                    startTime+adjmatrix[start][i]<= nodeGroup[i].closeTime ){
                    if (nodeGroup[i].happiness> nodeGroup[tmpSelect[0]].happiness){                    
                        tmpSelect.push_back(i);
                        tmpSelect.erase(tmpSelect.begin());
                    }
                }
                else if (startTime+adjmatrix[start][i] < nodeGroup[i].openTime){
                    waitTime=nodeGroup[i].openTime-(startTime+adjmatrix[start][i]);
                    if (travelTime-adjmatrix[start][i]-waitTime >= dist[i]){
                        if (nodeGroup[i].happiness> nodeGroup[tmpSelect[0]].happiness){                    
                            tmpSelect.push_back(i);
                            tmpSelect.erase(tmpSelect.begin());
                        }                                              

                    }
                }
            }
            else{
                //cout<<"empty"<<endl;
                findNext=true;           
                
                 if (startTime+adjmatrix[start][i]>= nodeGroup[i].openTime && 
                    startTime+adjmatrix[start][i]<= nodeGroup[i].closeTime ){
                        tmpSelect.push_back(i);
        
                }
                else if (startTime+adjmatrix[start][i] < nodeGroup[i].openTime){
                    waitTime=nodeGroup[i].openTime-(startTime+adjmatrix[start][i]);
                    if (travelTime-adjmatrix[start][i]-waitTime >= dist[i]){
                        tmpSelect.push_back(i);                                              
                    }
                }         
            }    
        }
    }

    if (tmpSelect.size()>0){
        bool stay=false;
        int next=tmpSelect[0];
        int remain=0;
        if (startTime+adjmatrix[start][next]>= nodeGroup[next].openTime && 
            startTime+adjmatrix[start][next]<= nodeGroup[next].closeTime ){
                stay=false;
                remain=travelTime-adjmatrix[start][next];
            }
        else{
            stay=true;
            remain=travelTime-nodeGroup[next].openTime+startTime;
        }
        stayStation.push_back(stay);       
        visited[next]=true;
        aPath.push_back(next);
        travelAns2(next, remain, nodeGroup, aPath, stayStation, startTime+travelTime-remain);
    }
    else if (findNext){  //have path no happy
        vector<int> tmpdist;
        for (int i=0; i<total_vertex; i++ ){
            if ((adjmatrix[start][i] !=INI_WEIGHT ) && (travelTime-adjmatrix[start][i] >= dist[i] ) && (visited[i]==false) ){
                if (tmpSelect.size()>0){
                    if (adjmatrix[start][i]< tmpdist[0]){                  
                        tmpSelect.push_back(i);
                        tmpSelect.erase(tmpSelect.begin());
                        tmpdist.push_back(adjmatrix[start][i]);
                        tmpdist.erase(tmpdist.begin());
                    }
                }
                else{
                    //cout<<"empty"<<endl;
                    tmpSelect.push_back(i);
                    tmpdist.push_back(adjmatrix[start][i]);
                }    
            }
        }
        int next=tmpSelect[0];
        int remain=travelTime-adjmatrix[start][next];
        visited[next]=true;
        aPath.push_back(next);
        stayStation.push_back(false);
        travelAns2(next, remain, nodeGroup, aPath, stayStation, startTime+travelTime-remain);    

    }
    else{
        gohomeAns2(start, travelTime,aPath, stayStation, nodeGroup, startTime);
    }
}

int GetIdxFromName(string aname, vector<nodeInfo> nodeGroup){
    int index=-1;
    nodeInfo tmpnode;
    for(int i=0; i<nodeGroup.size(); i++){
        tmpnode=nodeGroup[i];
        if (tmpnode.name==aname){
            index=tmpnode.graphIdx;
            break;
        }
    }
    if (index==-1){
        cout<<"Can not find node!"<<endl;
    }

    return index;
}


int main(int argc, char *argv[]){
//  read file
    string filePath;
    filePath.assign(argv[1]);
    filePath=filePath+"/tp.data";
    ifstream inputFile(filePath);
    if (!inputFile)
        cout<<"Can not open file!"<<endl;

    int nodeCount=0, edgeCount=0, travelTotalTime=0, travelStartTime=0; 
    string tmpstr,substr;

    int current=0, next=0;
    
        getline(inputFile, tmpstr);
         //for (int j=0; j<tmpstr.size(); j++)
            //cout<<tmpstr<<endl;
        current=0;
        next=0;
    

    for(int i=0 ;i<4;i++){
        next=tmpstr.find_first_of(" ",current);
        substr=tmpstr.substr(current, next-current);
        if (i==0)
            nodeCount=atoi(substr.c_str());
        else if(i==1)
            edgeCount=atoi(substr.c_str());
        else if(i==2)
            travelTotalTime=atoi(substr.c_str());
        else 
            travelStartTime=atoi(substr.c_str());
        
        current=next+1;
        
        
    }

    /*
    cout<<"done"<<endl;
    cout<<nodeCount<<endl;
    cout<<edgeCount<<endl;
    cout<<travelTotalTime<<endl;
    cout<<travelStartTime<<endl;*/

    tmpstr="";
    substr="";
    vector<nodeInfo> nodeGroup;
    nodeInfo tmpNode; 

    int happiness=0,openTime=0,closeTime=0;
    string nodeName="";
    for (int k=0; k<nodeCount; k++){
        getline(inputFile, tmpstr);
        //cout<<tmpstr<<endl;
        happiness=0,openTime=0,closeTime=0;
        nodeName="";
        current=0;
        for(int i=0 ;i<4;i++){
            next=tmpstr.find_first_of(" ",current);
            substr=tmpstr.substr(current, next-current);
            if (i==0)
                nodeName=substr;
            else if(i==1)
                happiness=atoi(substr.c_str());
            else if(i==2)
                openTime=atoi(substr.c_str());
            else 
                closeTime=atoi(substr.c_str());
            
            current=next+1;           
        } 
        tmpNode.name=nodeName;
        tmpNode.happiness=happiness;
        tmpNode.openTime=openTime;
        tmpNode.closeTime=closeTime;
        tmpNode.graphIdx=k;
        nodeGroup.push_back(tmpNode);
    }

 /*   for(int i=0; i<nodeGroup.size(); i++){
        tmpNode=nodeGroup[i];
        cout<<tmpNode.name<<" "<<tmpNode.happiness<<" "<<tmpNode.openTime<<" "<<tmpNode.closeTime<<endl;;
    }*/

    Graph travelMap(nodeCount);
    
    string node1,node2;
    int weight,idx1,idx2;
    for (int k=0; k<edgeCount; k++){
        tmpstr="";
        node1="";
        node2="";
        weight=0;
        idx1=-1;
        idx2=-1;
        getline(inputFile, tmpstr);
  
        happiness=0,openTime=0,closeTime=0;
        nodeName="";
        current=0;
        for(int i=0 ;i<3;i++){
            next=tmpstr.find_first_of(" ",current);
            substr=tmpstr.substr(current, next-current);
            if (i==0)
                node1=substr;
            else if(i==1)
                node2=substr;
            else 
                weight=atoi(substr.c_str());
           
            current=next+1;  
        }    

        idx1=GetIdxFromName(node1, nodeGroup);
        idx2=GetIdxFromName(node2, nodeGroup);

        //cout<<idx1<<" "<<idx2<<" "<<weight<<endl;

        travelMap.AddEdge(idx1,idx2, weight);


    }

    

    vector<int> tmpPath1, optPath1;
    int optHappy1=0, tmpHappy1=0;

    vector<int> tmpPath2, optPath2;
    int optHappy2=0, tmpHappy2=0;

    vector<bool> stayGroup,optStay;
    int waitTime=0;
    bool stay=false;
    int remainTime=0,startTime=0;

    for (int source=0;source<nodeCount;source++){
        waitTime=0;
        stay=false;
        remainTime=travelTotalTime;
        startTime=travelStartTime;
        vector<int>().swap(tmpPath1);
        vector<int>().swap(tmpPath2);
        vector<bool>().swap(stayGroup);
    

        travelMap.InitPredDist();
        travelMap.CalculateMinPath(source);
        travelMap.InitVisited();
        //cout<<"begin travel"<<endl;
        tmpPath1.push_back(source);
        travelMap.travel(source, travelTotalTime, nodeGroup, tmpPath1);

        travelMap.InitVisited();
        tmpPath2.push_back(source);
        if (travelStartTime<nodeGroup[source].openTime){
            waitTime=nodeGroup[source].openTime-travelStartTime;
            if (waitTime<=travelTotalTime){
                stay=true;
                startTime=startTime+waitTime;
                remainTime=remainTime-waitTime;
            }
            else{
                stay=false;
            }
        }
        else{
            stay=false;
        }        

        stayGroup.push_back(stay);
        travelMap.travelAns2(source, remainTime, nodeGroup, tmpPath2, stayGroup, startTime);

/*
        cout<<"stay size: "<<stayGroup.size()<<endl;
        for(int i=0; i<stayGroup.size(); i++){
            cout<<stayGroup[i]<<" ";
        }
        cout<<endl;     
    
        cout<<"ans1:"<<endl;
        for(int i=0; i<tmpPath1.size(); i++){
            cout<<tmpPath1[i]<<" ";
        }
        cout<<endl;    


        cout<<"ans2:"<<endl;
        for(int i=0; i<tmpPath2.size(); i++){
            cout<<tmpPath2[i]<<" ";
        }
        cout<<endl;       */ 
        

        tmpHappy1=travelMap.CalculateHappyAns1(nodeGroup,tmpPath1);

        tmpHappy2=travelMap.CalculateHappyAns2(nodeGroup, tmpPath2,stayGroup, startTime);

        //cout<<"happy2: "<<tmpHappy2<<endl;



        if (tmpHappy2>optHappy2){
            optHappy2=tmpHappy2;
            optPath2=tmpPath2;
            optStay=stayGroup;
        }
            
        
        if (tmpHappy1>optHappy1){
            optHappy1=tmpHappy1;
            optPath1=tmpPath1;
        }   

    }

      /*for(int i=0; i<optPath.size(); i++){
            cout<<optPath[i]<<" ";
        }
        cout<<"happy: "<<optHappy<<endl;*/

    vector<travelStep> stepGroup; 
    vector<travelStep> stepGroupAns2;

    int totalTimeCost1=0;
    travelStep tmpstep;
    totalTimeCost1=travelMap.CalculateTravelTimeAns1(nodeGroup,stepGroup, optPath1, travelStartTime);

    

    int totalTimeCostAns2=travelMap.CalculateTravelTimeAns2(nodeGroup,stepGroupAns2, optPath2, travelStartTime, optStay);
  /*  cout<<"time cost: "<<totalTimeCostAns2<<endl;
    for (int i=0; i<stepGroupAns2.size(); i++){
        tmpstep=stepGroupAns2[i];
        cout<<tmpstep.name<<" "<<tmpstep.enterTime<<" "<<tmpstep.leaveTime<<endl;
    } */  

 
    /*cout<<"time cost: "<<totalTimeCost<<endl;
    for (int i=0; i<stepGroup.size(); i++){
        tmpstep=stepGroup[i];
        cout<<tmpstep.name<<" "<<tmpstep.enterTime<<" "<<tmpstep.leaveTime<<endl;
    }*/

    string ans1path;
    ans1path.assign(argv[1]);
    ans1path=ans1path+"/ans1.txt";
    ofstream answerOutput (ans1path);

    if (answerOutput.is_open()){        
        answerOutput<<optHappy1<<" "<<totalTimeCost1<<endl;
        for (int i=0; i<stepGroup.size(); i++){
            tmpstep=stepGroup[i];
            answerOutput<<tmpstep.name<<" "<<tmpstep.enterTime<<" "<<tmpstep.leaveTime<<endl;
        }    
    }
    else cout<<"unable to open file!";   


    string ans2path;
    ans2path.assign(argv[1]);
    ans2path=ans2path+"/ans2.txt";
    ofstream ans2Output (ans2path);

    if (ans2Output.is_open()){
        ans2Output<<optHappy2<<" "<<totalTimeCostAns2<<endl;
        for (int i=0; i<stepGroupAns2.size(); i++){
            tmpstep=stepGroupAns2[i];
            ans2Output<<tmpstep.name<<" "<<tmpstep.enterTime<<" "<<tmpstep.leaveTime<<endl;
        }

    }
    else cout<<"unable to open file!";   




    return 0;
} 