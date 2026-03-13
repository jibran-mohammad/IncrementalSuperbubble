#include"superbubble.h"

    void SuperbubbleCalculation::computeSB(DFSCONTEXT *sdfs,GraphContext *g){
	    int t1,t;
	    std::stack<ppii> S;
	    for(t1=0;t1<=g->displayNumberVertices()-1;t1++){
	    	t=sdfs->dfsStrategy->topOrd[t1];
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

    
    void SuperbubbleCalculation::printSuperbubbles() {
    	std::cout << "superbubbles of graph are" << std::endl;
    	for (pii element : this->SB) {
            std::cout << "<"<<element.first <<" " << element.second <<">";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


    int SuperbubbleCalculation::numberOfSuperbubbles() {
    	int numElements = this->SB.size();
    	return numElements;
    }

    void SuperbubbleCalculation::insertE(DFSCONTEXT *sdfs,GraphContext *g,int u, int v) {
    	sdfs->updateEdge(g,u,v);
        this->SB.clear();
    }






