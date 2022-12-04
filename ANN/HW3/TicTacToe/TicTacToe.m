clc
clear

alpha = 0.1;
gamma = 1;
epsilon = 1;
p1_wins = 0;
p2_wins = 0;
draws = 0;

p1_wins_list = zeros([1,300]);
p2_wins_list= zeros([1,300]);
draws_list = zeros([1,300]);

Q1_table = cell(0);
Q2_table = cell(0);

for i = 1:3000
    
    if ~mod(i,100)
        epsilon = epsilon * 0.95;
        p1_wins_list(fix(i/100)) = p1_wins;
        p2_wins_list(fix(i/100)) = p2_wins;
        draws_list(fix(i/100)) = draws;
        p1_wins = 0;
        p2_wins = 0;
        draws = 0;

    end
    
    board = zeros(3);
    isend = 0;
    player = 1;

    Q1_table = ifseen(board, Q1_table);
    [action_board, action_coord_1] = choose_action(board, epsilon, player, Q1_table);
    board_prime = board + action_board;
    player = ~player;
    
    Q2_table = ifseen(board_prime, Q2_table);
    [action_board, action_coord_2] = choose_action(board_prime, epsilon, player, Q2_table);
    board_next = board_prime + action_board;
    
    Q1_table = ifseen(board_next, Q1_table);
    Q1_table = update_Q(board, board_next, action_coord_1, Q1_table, alpha, gamma, 0);
    player = ~player;
    board = board_next;
    while true
        board_prime_old = board_prime;
        [action_board, action_coord_1] = choose_action(board, epsilon, player, Q1_table);
        board_prime = board + action_board;
        isend = boardcheck(board_prime);
        if isend
            [R_1, R_2] = feed_reward(isend, player);
            Q1_table = fupdate_Q(board, action_coord_1, Q1_table, alpha, gamma, R_1);
            Q2_table = fupdate_Q(board_prime_old, action_coord_2, Q2_table, alpha, gamma, R_2);
            if R_1 == 1
                p1_wins = p1_wins + 1;
            elseif R_1 == 0
                draws = draws + 1;
            end
            break
        end
        Q2_table = ifseen(board_prime, Q2_table);
        Q2_table = update_Q(board_prime_old, board_prime, action_coord_2, Q2_table, alpha, gamma, 0);   
        player = ~player;
        
        [action_board, action_coord] = choose_action(board_prime, epsilon, player, Q2_table);
        board_next = board_prime + action_board;
        isend = boardcheck(board_next);
        if isend
            [R_1, R_2] = feed_reward(isend, player);
            Q1_table = fupdate_Q(board, action_coord_1, Q1_table, alpha, gamma, R_1);
            Q2_table = fupdate_Q(board_prime, action_coord_2, Q2_table, alpha, gamma, R_2);
            if R_2 == 1
                p2_wins = p2_wins + 1;
            elseif R_2 == 0
                draws = draws + 1;
            end
            break
        end
        Q1_table = ifseen(board_next, Q1_table);
        Q1_table = update_Q(board, board_next, action_coord_1, Q1_table, alpha, gamma, 0);
        player = ~player;
        board = board_next;
    end
end

figure
plot(p1_wins_list)
hold on
plot(p2_wins_list)
plot(draws_list)
legend("P1 wins", "P2 wins", "Draws")


%Functions:

%Reward
function [R_1, R_2] = feed_reward(isend, player)
    R_1 = 0;
    R_2 = 0;
    if isend == 1 && player == 1
        R_1 = 1;
        R_2 = -1;
    elseif isend == 1 && player == 2
        R_1 = -1;
        R_2 = 1;
    elseif isend == 0.5
        R_1 = 0;
        R_2 = 0;
    end
end

function Q_table = ifseen(board, Q_table)
    if isempty(Q_table)

        Q_table{1,1} = board;
        board(board~=0)=NaN;
        Q_table{2,1} = board;
    else
        index = find(cellfun(@(x) isequal(x, board), Q_table(1,:)), 1);
        if isempty(index)
            
            Q_table{1,end+1} = board;
            board(board~=0)=NaN;
            Q_table{2,end} = board;
        end
    end  
end



function [action_board, action_coord] = choose_action(board, epsi, which_player, Q_table)
    action_board = zeros(3);
    if binornd(1, epsi) == 1
        [row, col] = find(board==0);
        randchoice = randi(length(row));
        action_coord = [row(randchoice); col(randchoice)];
    else

        index = find(cellfun(@(x) isequal(x, board), Q_table(1,:)), 1);
        max_qvalue = max(Q_table{2, index}(:));

        [x,y] = find(Q_table{2, index} == max_qvalue);
        if length(x) > 1
            randchoice = randi(length(x));
            action_coord = [x(randchoice); y(randchoice)];
        else
            action_coord = [x;y];
        end
    end
    if which_player == 1
        action_board(action_coord(1),action_coord(2)) = 1;
    else
        action_board(action_coord(1),action_coord(2)) = -1;
    end
end

%% Q update when game continues
function Q_table = update_Q(board_old, board_new, action_coord, Q_table, alpha, gamma, reward)
    index_old = find(cellfun(@(x) isequal(x, board_old), Q_table(1,:)), 1);
    index_new = find(cellfun(@(x) isequal(x, board_new), Q_table(1,:)), 1);
    max_Q = max(Q_table{2, index_new}(:));
    Q_table{2, index_old}(action_coord(1), action_coord(2)) = Q_table{2, index_old}(action_coord(1), action_coord(2)) + ...
                                                                alpha * (reward + gamma * max_Q - Q_table{2, index_old}(action_coord(1), action_coord(2)));
end

%% Q update when game over
function Q_table = fupdate_Q(board_old, action_coord, Q_table, alpha, gamma, reward)
    index_old = find(cellfun(@(x) isequal(x, board_old), Q_table(1,:)), 1);
    Q_table{2, index_old}(action_coord(1), action_coord(2)) = Q_table{2, index_old}(action_coord(1), action_coord(2)) + ...
                                                                alpha * (reward - Q_table{2, index_old}(action_coord(1), action_coord(2)));
end

%% Check if game over
function isOver = boardcheck(board)


    persistent p1_winning;
    persistent p2_winning;
    p1_winning = [1,1,1];
    p2_winning = [-1,-1,-1];

    if ~ismember(0,board)
        isOver = 0.5;
        return
    end

    
    v = ismember(board, p1_winning,'rows');
    r = find(v,1);
    if r>0
        isOver = 1;
        return
    end
    v = ismember(board',p1_winning, "rows");
    r = find(v,1);
    if r>0
        isOver = 1;
        return
    end
    v = ismember(board, p2_winning,'rows');
    r = find(v,1);
    if r>0
        isOver = 1;
        return
    end
    v = ismember(board',p2_winning, "rows");
    r = find(v,1);
    if r>0
        isOver = 1;
        return
    end
    
    
    
    elem = board(1,1);
    if elem ~=0
        for i=2:3
            if board(i,i)==elem
                if i==3
                    isOver = 1;
                    return
                end
            else
                break
            end
        end
    end
    
    elem = board(1,3);
    if elem ~= 0
        for i=2:3
            if board(i,4-i)==elem
                if i==3
                    isOver = 1;
                    return
                end
            else
                break
            end
        end
    end
    isOver = 0;
end