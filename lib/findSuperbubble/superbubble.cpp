#include"superbubble.h"

    void SuperbubbleCalculation::computeSB(DFSCONTEXT *sdfs,GraphContext *g){
	    int t1,t;
	    std::stack<ppii> S;
	    for(t1=0;t1<=g->displayNumberVertices()-1;t1++){
	    	t=sdfs->dfsStrategy->topOrd[t1];
	    	//std::cout << t << std::endl;
	    	/*
	    	if(sdfs->dfsStrategy->minI[t] == INT_MAX) {
	    		while(!S.empty()){
	    			S.pop();
	    		}
	    	}
	    	*/
	    	// if minI[t] goes beyond what is on top of the stack, then pop elements until the elemement is found
	    	// which does where minI[t] is within the bounds.
		    while(!S.empty() && sdfs->dfsStrategy->reverseTopologicalOrder[S.top().first] > sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]]){
			    ppii s1 = S.top();
			    S.pop();
	            if(!S.empty()) {
			        ppii s2 = S.top();
			        S.pop();
			        if(s1.second.first == INT_MAX) {
			        	S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(s1.second.first, sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			        else {
			        S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first], sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			    }

		    }
			
		    if(!S.empty() && S.top().second.second == t){
			    //assert(MIN(S.top().second.first,sdfs->dfsStrategy->minI[t])==S.top().first);
			    this->SB.push_back(pii(S.top().first,t));
			    //printf("<%d,%d>\n", S.top().first,t);
			    ppii s1 = S.top();
			    S.pop();
			    if(!S.empty()){
				    ppii s2 = S.top();
				    S.pop();
			        if(s1.second.first == INT_MAX) {
			        	S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(s1.second.first, sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			        else {
			            S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first], sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
				       
			       
			    }

			    
		    }
		
		    if(!S.empty()){
			    ppii s1 = S.top();
			    S.pop();
		        if(s1.second.first == INT_MAX) {
		        	s1.second.first = sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]],s1.second.first)];
		        }
		        else {
			        s1.second.first = sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]],sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first])];
			    }
			    S.push(s1);
		    }
		
		    S.push(ppii(t,pii(INT_MAX,sdfs->dfsStrategy->maxO[t])));
	    }	

    }

    
    // SB holds internal vertex ids; they are converted to input ids only for printing.
    void SuperbubbleCalculation::printSuperbubbles(GraphContext *g) {
    	std::cout << "superbubbles of graph are" << std::endl;
    	for (pii element : this->SB) {
            std::cout << "<"<<g->toOriginal(element.first) <<" " << g->toOriginal(element.second) <<">";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


    // One superbubble per line, "<entrance> <exit>", in input ids.
    void SuperbubbleCalculation::writeSuperbubbles(GraphContext *g, std::ostream &out) {
        for (pii element : this->SB) {
            out << g->toOriginal(element.first) << " " << g->toOriginal(element.second) << "\n";
        }
    }


    int SuperbubbleCalculation::numberOfSuperbubbles() {
    	int numElements = this->SB.size();
    	return numElements;
    }

    void SuperbubbleCalculation::insertE(DFSCONTEXT *sdfs,GraphContext *g,int u, int v) {
    	sdfs->updateEdge(g,u,v);
        this->SB.clear();
    }
 
    // a batch of edges (input ids): the DFS strategy updates the order, minI
    // and maxO once for all of them
    void SuperbubbleCalculation::insertBatch(DFSCONTEXT *sdfs, GraphContext *g, const std::vector<std::pair<int,int>> &edges) {
	    sdfs->updateEdgeBatch(g, edges);
	    this->SB.clear();
	}






