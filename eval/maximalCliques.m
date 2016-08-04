function complete_subgraphs = maximalCliques( X )
%MAXIMALCLIQUES finds all the maximal complete sub-graphs in a graph
%   The graph passed must be an upper rectangular square matrix. Each row
%   of a matrix denotes the presence of an edge with 1, and an absence of
%   it with 0. The row and col no. of an edge denotes the connecting nodes.
%   Given this matrix, this function finds all the maximal complete 
%   sub-graph (a set of nodes amongst all the nodes which form a complete 
%   sub-graph i.e. every node connects to every other) also known as 
%   cliques. A maximal graph is returned since every complete sub-graph 
%   will also have smaller complete sub-graphs inside itself. NOTE: this 
%   function would not return single node sub-graphs, although every
%   isolated node, in concept, also forms a complete sub-graph
%   The function returns all the sub-graphs in a cell-array, where each
%   row denotes a new sub-graph

%TEST CASES

% A = [-1  1  1  0  0  0
%      -1 -1  1  0  0  0
%      -1 -1 -1  1  0  1 
%      -1 -1 -1 -1  0  1
%      -1 -1 -1 -1 -1  0
%      -1 -1 -1 -1 -1 -1 ];

% B = [-1  1  1  0  1  1
%      -1 -1  1  1  1  1
%      -1 -1 -1  1  1  1 
%      -1 -1 -1 -1  1  0
%      -1 -1 -1 -1 -1  1
%      -1 -1 -1 -1 -1 -1 ];

    [m n] = size(X);
    
    assert( m == n, 'The matrix should be square, cause each side denotes the nodes in the same graph' );

    %the graph 
    complete_subgraphs = {};
    
    %this will keep track of all the nodes that need to be looked at
    remaining_nodes = 1:n;
    
    %discard the lower triangular matrix, just in case it has values
    X = X - tril(X);
    
    %main loop
    while ~isempty(remaining_nodes)
        
        %choose just the first node, in the possible nodes
        source_node = remaining_nodes(1);
        
        %remove from remaining_nodes
        remaining_nodes = setdiff(remaining_nodes, source_node);
        
        %make new set for the node
        new_complete_subgraph = source_node;
        
        %find first node with whom the source node has a edge
        connecting_nodes = returnConnectingNodes( X, source_node );
        
        %only if the node is not isolated
        if ~isempty(connecting_nodes)
            %add first connecting node to the the subgraph
            new_complete_subgraph = [ source_node connecting_nodes(1) ];
            
            %remove newly added node from remaining nodes
            remaining_nodes = setdiff(remaining_nodes, connecting_nodes(1));
            
            %make a list of remaining nodes which needs to be still to be explored
            other_psb_nodes = setdiff( connecting_nodes, new_complete_subgraph );

            %form a maximal sub-graph
            for i = other_psb_nodes
                %check if the node can be added to the current complete subgraph
                if checkNodeForMaximalSubgraph( X, new_complete_subgraph, i ) == 1
                    
                    %add node to the complete subgraph
                    new_complete_subgraph = [ new_complete_subgraph i ];
                    
                    %remove newly added node from remaining nodes
                    remaining_nodes = setdiff(remaining_nodes, i);
                end
            end
        
            %add subgraph to other subgraphs
            complete_subgraphs = [ complete_subgraphs; sort(new_complete_subgraph) ];
        end
        
    end
    
end


function can_be_added = checkNodeForMaximalSubgraph( X, complete_subgraph, test_node )
%Returns if the given test_node has edges to all the nodes in the
%complete_subgraph

    %find the connecting nodes of the test node
    test_node_edges = returnConnectingNodes( X, test_node );
    
    %see if all the nodes already present in the subgraph, are also present
    %in the connecting nodes of the test node
    can_be_added = all( ismember( complete_subgraph, test_node_edges ) );
end


function connecting_nodes = returnConnectingNodes( X, node_no )
%Returns the nodes to which the node_no have edges

    row_psbs = find( X(node_no,:) );        %search row for psbs
    col_psbs = find( X(:,node_no) );        %search col for psbs
    
    %combine all the edges and return
    connecting_nodes = union( row_psbs, col_psbs );
end