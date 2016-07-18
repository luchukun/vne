//min_cost flow by network simplex algorithm
#ifndef MIN_COST_FLOW_H
#define MIN_COST_FLOW_H
#include <vector>
#include <limits>
#include <algorithm>
class network_simplex{
      
	public:
	enum PivotRule {
	  /// The \e First \e Eligible pivot rule.
	  /// The next eligible arc is selected in a wraparound fashion
	  /// in every iteration.
	  FIRST_ELIGIBLE,

	  /// The \e Best \e Eligible pivot rule.
	  /// The best eligible arc is selected in every iteration.
	  BEST_ELIGIBLE,

	  /// The \e Candidate \e List pivot rule.
	  /// In a major iteration a candidate list is built from eligible arcs
	  /// in a wraparound fashion and in the following minor iterations
	  /// the best eligible arc is selected from this list.
	  CANDIDATE_LIST,
	};
	// State constants for arcs
	enum ArcState {
	  STATE_UPPER = -1,
	  STATE_TREE  =  0,
	  STATE_LOWER =  1
	};

	// Direction constants for tree arcs
	enum ArcDirection {
	  DIR_DOWN = -1,
	  DIR_UP   =  1
	};

// for each edge
	vector<int> _source;
	vector<int> _target;
	vector<float> _capacity;
	vector<float> _flow;
	vector<float> _cost;
	vector<int> _state;
// for each node
	vector<float> _potential;
	vector<int> _parent;
	vector<int> _depth;
	vector<int> _thread;
	int _root;
	int _next_arc;
	int _enter_arc;
	int apex;
	int _search_arc_num;

	// Implementation of the First Eligible pivot rule
	bool FirstEligible(int rule ) 
	{	
		float c; //reduced cost
		for (int e = _next_arc; e != _search_arc_num; ++e) {
		  c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
		  if (c < 0) {
			_enter_arc = e;
			_next_arc = e + 1;
			return true;
		  }
		}
		for (int e = 0; e != _next_arc; ++e) {
		  c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
		  if (c < 0) {
			_enter_arc = e;
			_next_arc = e + 1;
			return true;
		  }
		}
		return false;
	}
	// Implementation of the Best Eligible pivot rule	        
	bool BestEligible() {
		   float c, min = 0; //reduced cost
			for (int e = 0; e != _search_arc_num; ++e) {
			  c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
			  if (c < min) {
				min = c;
				_enter_arc = e;
			  }
			}
			return min < 0;
		  }


	// Implementation of the Candidate List pivot rule
	bool Candidate List() {
			Cost min, c;
			int e;
			if (_curr_length > 0 && _minor_count < _minor_limit) {
			  // Minor iteration: select the best eligible arc from the
			  // current candidate list
			  ++_minor_count;
			  min = 0;
			  for (int i = 0; i < _curr_length; ++i) {
				e = _candidates[i];
				c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
				if (c < min) {
				  min = c;
				  _enter_arc = e;
				}
				else if (c >= 0) {
				  _candidates[i--] = _candidates[--_curr_length];
				}
			  }
			  if (min < 0) return true;
			}

			// Major iteration: build a new candidate list
			min = 0;
			_curr_length = 0;
			for (e = _next_arc; e != _search_arc_num; ++e) {
			  c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
			  if (c < 0) {
				_candidates[_curr_length++] = e;
				if (c < min) {
				  min = c;
				  _enter_arc = e;
				}
				if (_curr_length == _list_length) goto search_end;
			  }
			}
			for (e = 0; e != _next_arc; ++e) {
			  c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
			  if (c < 0) {
				_candidates[_curr_length++] = e;
				if (c < min) {
				  min = c;
				  _enter_arc = e;
				}
				if (_curr_length == _list_length) goto search_end;
			  }
			}
			if (_curr_length == 0) return false;

		  search_end:
			_minor_count = 1;
			_next_arc = e;
			return true;
		  }

	// Find the apex node
	int findJoinNode() {
	  int u = _source[enter_arc];
	  int v = _target[enter_arc];

	  while (u != v) {
		if (_depth[u] == _depth[v]) {
		  u = _parent[u];
		  v = _parent[v];
		}
		else if (_depth[u] > _depth[v]) u = _parent[u];
		else v = _parent[v];
	  }
	  return u;
	}

	// Find the leaving arc of the cycle and returns true if the
	// leaving arc is not the same as the entering arc
	bool findLeavingArc(int enter_arc,int apex,float& delta) {
	  // Initialize first and second nodes according to the direction
	  // of the cycle
	  int first, second;
	  if (_state[enter_arc] == STATE_LOWER) {
		first  = _source[enter_arc];
		second = _target[enter_arc];
	  } else {
		first  = _target[enter_arc];
		second = _source[enter_arc];
	  }
	  delta = _capacity[enter_arc];
	  bool result = false;
	  float d;
	  int e; //edge

	  // Searching the cycle along the path form the first node to the
	  // root node
	  for (int u = first; u != apex; u = _parent[u]) {
		e = _pred_edge[u];
		d = _forward[u] ? _flow[e] : _capacity[e] - _flow[e];
		if (d < delta) {
		  delta = d;
		  u_out = u; //leaving arc?
		  u_in = first;
		  v_in = second;
		  result = true;
		}
	  }
	  // Searching the cycle along the path form the second node to the
	  // root node
	  for (Node u = second; u != apex; u = _parent[u]) {
		e = _pred_edge[u];
		d = _forward[u] ? _capacity[e] - _flow[e] : _flow[e];
		if (d <= delta) {
		  delta = d;
		  u_out = u;
		  u_in = second;
		  v_in = first;
		  result = true;
		}
	  }
	  return result;
	}

	/// Changes \c flow and \c state edge maps.
	void changeFlows(bool change,float delta) {
	  // Augmenting along the cycle
	  if (delta > 0) {
		float val = _state[_enter_arc] * delta;
		_flow[_enter_arc] += val;
		for (Node u = _graph.source(_enter_arc); u != apex; u = _parent[u]) {
		  _flow[_pred_edge[u]] += _forward[u] ? -val : val;
		}
		for (Node u = _graph.target(_enter_arc); u != apex; u = _parent[u]) {
		  _flow[_pred_edge[u]] += _forward[u] ? val : -val;
		}
	  }
	  // Updating the state of the entering and leaving edges
	  if (change) {
		_state[_enter_arc] = STATE_TREE;
		_state[_pred_edge[u_out]] =
		  (_flow[_pred_edge[u_out]] == 0) ? STATE_LOWER : STATE_UPPER;
	  } else {
		_state[_enter_arc] = -_state[_enter_arc];
	  }
	}

	/// Updates \c thread and \c parent node maps.
	void updateThreadParent() {
	  Node u;
	  v_out = _parent[u_out];

	  // Handling the case when apex and v_out coincide
	  bool par_first = false;
	  if (apex == v_out) {
		for (u = apex; u != u_in && u != v_in; u = _thread[u]) ;
		if (u == v_in) {
		  par_first = true;
		  while (_thread[u] != u_out) u = _thread[u];
		  first = u;
		}
	  }

	  // Finding the last successor of u_in (u) and the node after it
	  // (right) according to the thread index
	  for (u = u_in; _depth[_thread[u]] > _depth[u_in]; u = _thread[u]) ;
	  right = _thread[u];
	  if (_thread[v_in] == u_out) {
		for (last = u; _depth[last] > _depth[u_out]; last = _thread[last]) ;
		if (last == u_out) last = _thread[last];
	  }
	  else last = _thread[v_in];

	  // Updating stem nodes
	  _thread[v_in] = stem = u_in;
	  par_stem = v_in;
	  while (stem != u_out) {
		_thread[u] = new_stem = _parent[stem];

		// Finding the node just before the stem node (u) according to
		// the original thread index
		for (u = new_stem; _thread[u] != stem; u = _thread[u]) ;
		_thread[u] = right;

		// Changing the parent node of stem and shifting stem and
		// par_stem nodes
		_parent[stem] = par_stem;
		par_stem = stem;
		stem = new_stem;

		// Finding the last successor of stem (u) and the node after it
		// (right) according to the thread index
		for (u = stem; _depth[_thread[u]] > _depth[stem]; u = _thread[u]) ;
		right = _thread[u];
	  }
	  _parent[u_out] = par_stem;
	  _thread[u] = last;

	  if (apex == v_out && par_first) {
		if (first != v_in) _thread[first] = right;
	  } else {
		for (u = v_out; _thread[u] != u_out; u = _thread[u]) ;
		_thread[u] = right;
	  }
	}

	/// Updates \c pred_edge and \c forward node maps.
	void updatePredEdge() {
	  Node u = u_out, v;
	  while (u != u_in) {
		v = _parent[u];
		_pred_edge[u] = _pred_edge[v];
		_forward[u] = !_forward[v];
		u = v;
	  }
	  _pred_edge[u_in] = _enter_arc;
	  _forward[u_in] = (u_in == _graph.source(_enter_arc));
	}

	/// Updates \c depth and \c potential node maps.
	void updateDepthPotential() {
	  _depth[u_in] = _depth[v_in] + 1;
	  _potential[u_in] = _forward[u_in] ?
		_potential[v_in] - _cost[_pred_edge[u_in]] :
		_potential[v_in] + _cost[_pred_edge[u_in]];

	  Node u = _thread[u_in], v;
	  while (true) {
		v = _parent[u];
		if (v == INVALID) break;
		_depth[u] = _depth[v] + 1;
		_potential[u] = _forward[u] ?
		  _potential[v] - _cost[_pred_edge[u]] :
		  _potential[v] + _cost[_pred_edge[u]];
		if (_depth[u] <= _depth[v_in]) break;
		u = _thread[u];
	  }
	}

	/// Executes the algorithm.
	/// \brief Extends the underlaying graph and initializes all the
	/// node and edge maps.
	bool init(Graph &_graph,int e,float* hostBw) 
	{
	int s,d;
	float TS[Nv]={0},TD[Nv]={0};
	float TS_max[Nv]={0},TD_max[Nv]={0};
	int source_idx[Nv]={0},dest_idx[Nv]={0};
	float T_left=0,T_add=0,T_max=0;
	float sum_t=0; // sum of traffic on link e(i,j)
	if(_graph.n_pair[e]==0)
		return 0;
	int num_arc=0,num_vertex=1;
	//---- Initialize max concurrent traffic T_max[s] for each server s
	for(int p=0;p<_graph.n_pair[e];p++)
	{
		s = _graph.pass[e][p].s; d = _graph.pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			//从s出发的总流量限制
			TS_max[s] = hostBw[s];
			//到达d的总流量限制
			TD_max[d] = hostBw[d];
			num_arc++;
		}
	}
	for (int i=0;i<Nv;i++)
	{	if (TS_max[i]>0)
		{	source_idx[i]=num_vertex;
			num_vertex++;
			num_arc++;
		}
	}
	for (int i=0;i<Nv;i++)
	{	if (TD_max[i]>0)
		{	dest_idx[i]=num_vertex;
			num_vertex++;
			num_arc++;				
		}
	}
			
	// for each edge
	_source.resize(num_arc);
	_target.resize(num_arc);
	_capacity.resize(num_arc);
	_flow.resize(num_arc);
	_cost.resize(num_arc);
	_state.resize(num_arc);
// for each node
	_potential.resize(num_vertex);
	_parent.resize(num_vertex,-1);
	_depth.resize(num_vertex,-1);
	_thread.resize(num_vertex,-1);
	
	int r=0; //index of arc
	for(int p=0,e=0;p<_graph.n_pair[e];p++)
	{
		s = _graph.pass[e][p].s; d = _graph.pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			// compute a feasible flow by greedy method
			T_left = min(TS_max[s]-TS[s], TD_max[d]-TD[d]);
			T_add =max((float)0, min( min(hostBw[s],hostBw[d]), T_left));
			sum_t += T_add*_graph.f[s][d][e]; //max traffic
			TS[s] += T_add;	
			TD[d] += T_add;

			// variable initiation
			_source[r]=source_idx[s];_target[r]=dest_idx[d];
			_capactiy[r]=min(hostBw[s],hostBw[d]);
			_cost[r]=_graph.f[s][d][e];
			_flow[r]=T_add;			
			r++;
		}
	}
	_root=0;
	for (int i=0;p<Nv;p++)
	{	
		if (TS_max[i]>0)
		{	_source[r]=_root; _target[r]=source_idx[i]; 
			_flow[r]=TS[i];
			_capacity[r]=TS_max[i];
			_cost[r]=0;
			r++;
		}
		if (TD_max[i]>0)
		{	_source[r]=dest_idx[i]; _target[r]=_root; 
			_flow[r]=TD[i];
			_capacity[r]=TD_max[i];
			_cost[r]=0;
			r++;
		}
	}
	vector<bool>connected; connected.resize(num_vertex);

	for (r=0;r<num_arc;r++)
	{
		if(_flow[r]==0)
			_state[r]=STATE_LOWER;
		else if(_flow[r]==_capacity[r])
			_state[r]=STATE_UPPER;
		else
		{
			_state[r]=STATE_TREE;
			connected[_source[r]]=true;
			connected[_target[r]]=true;
		}
	}

	for (r=0;r<num_arc;r++)
	{

	}	  

	  return true;
	}

	bool start() {

	  // Executing the network simplex algorithm
	  while (findEnteringEdge()) {
		apex = findJoinNode();
		bool change = findLeavingEdge();
		changeFlows(change);
		if (change) {
		  updateThreadParent();
		  updatePredEdge();
		  updateDepthPotential();
		}
	  }

	  // Checking if the flow amount equals zero on all the artificial
	  // edges
	  for (InEdgeIt e(_graph, _root); e != INVALID; ++e)
		if (_flow[e] > 0) return false;
	  for (OutEdgeIt e(_graph, _root); e != INVALID; ++e)
		if (_flow[e] > 0) return false;

	  // Copying flow values to _flow_result
	  if (_lower) {
		for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e)
		  (*_flow_result)[e] = (*_lower)[e] + _flow[_edge_ref[e]];
	  } else {
		for (typename Graph::EdgeIt e(_graph_ref); e != INVALID; ++e)
		  (*_flow_result)[e] = _flow[_edge_ref[e]];
	  }
	  // Copying potential values to _potential_result
	  for (typename Graph::NodeIt n(_graph_ref); n != INVALID; ++n)
		(*_potential_result)[n] = _potential[_node_ref[n]];

	  return true;
	}


}
#endif