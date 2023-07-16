#ifndef _WIN32
#pragma GCC optimize("O3, inline,omit-frame-pointer")
#endif
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <map>
#include <list>
#include <float.h>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

constexpr int N{2}; //Number of players
constexpr int expertiseCoefficient{10}; //Expertise Weight

enum Location{SAMPLES=0,DIAGNOSIS=1,MOLECULES=2,LABORATORY=3,START=4};
enum ActionType{GOTO=0, CONNECT=1, WAIT=2};

const vector<string> locationToStr{"SAMPLES","DIAGNOSIS","MOLECULES","LABORATORY","START_POS"};
const map<string,Location> StrToLocation{{"SAMPLES",SAMPLES},{"DIAGNOSIS",DIAGNOSIS},{"MOLECULES",MOLECULES},{"LABORATORY",LABORATORY}};
const array<Location,4> intToLocation{SAMPLES,DIAGNOSIS,MOLECULES,LABORATORY};
const vector<string> typeToStr{"A","B","C","D","E"};
const map<string,int> TypeToInt{{"A",0},{"B",1},{"C",2},{"D",3},{"E",4}};
//Distances from everywhere to everywhere
constexpr array<array<int,4>,5> Distances{array<int,4>{0,3,3,3},array<int,4>{3,0,3,4},array<int,4>{3,3,0,3},array<int,4>{3,4,3,0},array<int,4>{2,2,2,2}};

struct Timeout {
    high_resolution_clock::time_point startTime;
    double maxTimeSeconds;

    Timeout(double maxTimeSeconds) {
        this->startTime = high_resolution_clock::now();
        this->maxTimeSeconds = maxTimeSeconds;
    }

    bool isElapsed()
    {
        duration<double> time_span = duration_cast<duration<double>>(high_resolution_clock::now() - startTime);
        return time_span.count() >= maxTimeSeconds;
    }
};

struct sample_template{
    array<int,5> Cost;
    int score,exp;
};

struct genetic_sample: sample_template{
    bool diagnosed;
    int id,rank,owner;
    inline void operator=(const sample_template &a)noexcept{
        Cost=a.Cost;
        score=a.score;
        exp=a.exp;
    }
    inline genetic_sample(const sample_template &a)noexcept{
        *this=a;
    }
};

struct Robot 
{
    Location location;
    int eta;
    int score;
    array<int,5> Molecules,Expertise;
    vector<genetic_sample> samples;
    int _id;

    int totalExpertise()
    {
        return Expertise[0] + Expertise[1] + Expertise[2] + Expertise[3] + Expertise[4];
    }

    int totalMolecules()
    {
        return Molecules[0] + Molecules[1] + Molecules[2] + Molecules[3] + Molecules[4];
    }

    bool hasSpace()
    {
        return totalMolecules() < 10;
    }
};

struct Project{
    array<int,5> Target;
};

struct Action
{
    ActionType type;
    int id;
};

struct variation{
    double score;
    vector<array<Action,2>> Moves;
};

struct State
{
    int sampleCount;
    int turn = 0;
    array<Robot, 2> players;
    array<int, 5> Available;
    
    vector<Project> projects;
    vector<genetic_sample> Samples;

    array<list<sample_template>, 3> SamplePool;

    //void update(const array<Action,N> &M);
};

bool Completed_Project(const Robot &p,const Project &proj){
    for(int m=0;m<5;++m){
        if(p.Expertise[m]<proj.Target[m]){
            return false;
        }
    }
    return true;
};

bool ReadyToProduce(const Robot &p,const genetic_sample &s){
    if(!s.diagnosed){
        return false;
    }
    for(int i=0;i<5;++i){
        if(p.Molecules[i]<s.Cost[i]-p.Expertise[i]){
            return false;
        }
    }
    return true;
}

bool MoleculesAvailable(const Robot &p, const State &state,const genetic_sample &s) {
    for(int i=0;i< 5;i++) {
        if ( s.Cost[i] - p.Expertise[i] <= p.Molecules[i] + state.Available[i] ) {
            return false;
        }
    }
    return true;
}

int MissingMolecules(const Robot &p, const genetic_sample &s){
    int missingMolecules = 0;
    for(int i=0;i<5;i++) {
        int missing = s.Cost[i] - p.Expertise[i] - p.Molecules[i];
        if ( missing > 0) {
            missingMolecules = missingMolecules + missing;
        }
    }
    return missingMolecules;
}

bool EnoughSpaceToTake(const Robot &p, const genetic_sample &s) {
    return 10 -(p.Molecules[0] + p.Molecules[1] + p.Molecules[2] + p.Molecules[3] + p.Molecules[4]) > MissingMolecules(p, s);
}

double sampleScore(const Robot &p, const State &state,const genetic_sample &s){
    bool moleculesAvailable = MoleculesAvailable(p,state,s);
    return (moleculesAvailable && EnoughSpaceToTake(p,s) ? 100 : moleculesAvailable ? 1 : 0)-MissingMolecules(p, s)*1e-3;
}

// I need to test this
bool CouldProduce(const Robot &p,const State &state, const genetic_sample &sample){
    if(!sample.diagnosed){
        return false;
    }
    int totalNeeded = 0;
    for(int i=0;i<5;++i){
        
        if(state.Available[i] < sample.Cost[i]-p.Expertise[i]){
            return false;
        }
        if ( sample.Cost[i] - p.Expertise[i] > 0) {
            totalNeeded = totalNeeded + sample.Cost[i] - p.Expertise[i];
        }
    }

    int totalStorage = p.Molecules[0] + p.Molecules[1] + p.Molecules[2] + p.Molecules[3] + p.Molecules[4];
    if ( 10 - totalStorage < totalNeeded ) {
        return false;
    }
    return true;
}

void Simulate(State &S,const array<Action,N> &M){
    const State S_before=S;
    
    for(int i=0;i<2;++i){
        Robot& p{S.players[i]};
        const Action& mv{M[i]};
        if(p.eta==0){//Ignore actions of moving players
            if(mv.type==GOTO){
                p.eta=Distances[p.location][mv.id];
                p.location=intToLocation[mv.id];
            }
            else{//Connect
                if(p.location  ==SAMPLES){//Take undiagnosed sample
                    const int& rank{mv.id};
                    S.SamplePool[rank-1].push_back(S.SamplePool[rank-1].front());//Make a copy of the taken sample at the back of the list
                    genetic_sample s=S.SamplePool[rank-1].front();
                    S.SamplePool[rank-1].pop_front();
                    s.id=S.sampleCount++;
                    s.rank=rank;
                    s.diagnosed=false;
                    s.owner=i;
                    p.samples.push_back(s);
                }
                else if(p.location==MOLECULES){//Take molecule
                    if(S_before.Available[mv.id]>0){
                        --S.Available[mv.id];
                        ++p.Molecules[mv.id];
                    }
                }
                else if(p.location==LABORATORY && find_if(p.samples.begin(),p.samples.end(),[&](const genetic_sample &s){return s.id==mv.id;})!=p.samples.end()){
                    const auto s=find_if(p.samples.begin(),p.samples.end(),[&](const genetic_sample &s){return s.id==mv.id;});
                    if(ReadyToProduce(p,*s)){
                        for(int m=0;m<5;++m){
                            const int spent{max(0,s->Cost[m]-p.Expertise[m])};
                            p.Molecules[m]-=spent;
                            S.Available[m]+=spent;
                        }
                        p.score+=s->score;//Increase score
                        ++p.Expertise[s->exp];//Gain expertise
                        p.samples.erase(s);
                    }
                    else{
                        cerr << i << " tried to produce something he can't" << endl;
                    }
                }
                else if(p.location==DIAGNOSIS){
                    const auto player_s{find_if(p.samples.begin(),p.samples.end(),[&](const genetic_sample &s){return s.id==mv.id;})};
                    const auto diag_s{find_if(S.Samples.begin(),S.Samples.end(),[&](const genetic_sample &s){return s.id==mv.id;})};
                    if(player_s!=p.samples.end()){
                        if(player_s->diagnosed){
                            S.Samples.push_back(*player_s);
                            p.samples.erase(player_s);
                        }
                        else{
                            player_s->diagnosed=true;
                        }
                    }
                    else if(diag_s!=S.Samples.end()){
                        const Action mv2{M[(i+1)%2]};
                        const Robot& p2{S.players[(i+1)%2]};
                        const bool other_hasnt_requested{mv.id!=mv2.id || mv.type!=mv2.type || p2.location!=DIAGNOSIS || p2.eta>0 || p2.samples.size()==3};
                        if(p.samples.size()<3 && (diag_s->owner==i || other_hasnt_requested) ){
                            //cerr << "Player " << i << " got sample " << mv.id << endl;
                            p.samples.push_back(*diag_s);
                            S.Samples.erase(diag_s);
                        }
                    }
                }
            }
        }
    }
    for(int i=0;i<2;++i){//Decrease eta of both players
        Robot& p{S.players[i]};
        p.eta=max(0,p.eta-1);
    }
    for(auto it=S.projects.begin();it!=S.projects.end();){
        bool completed{false};
        for(int i=0;i<2;++i){
            if(Completed_Project(S.players[i],*it)){
                //cerr << "Player " << i << " completed a project" << endl;
                S.players[i].score+=50;
                completed=true;
            }
        }
        if(completed){
            it=S.projects.erase(it);
        }
        else{
            ++it;
        }
    }
}

double evalSample(State state, Robot player, genetic_sample sample) {
    int minScore = 0;
    if ( sample.rank == 1 ) {
        minScore = 1;
    } else if ( sample.rank == 2) {
        minScore = 10;
    } else {
        minScore = 30;
    }

    if( !sample.diagnosed ) {
        return 0.15 * (minScore + expertiseCoefficient);
    } 

    //TODO 0.175*(MinScore[rank]+ExpertiseCoeff) for a diagnosed unknown sample How to identify the simulated sample

    if( CouldProduce(player, state, sample) ){
        return 0.5*(sample.score + expertiseCoefficient);
    } else {
        return 0.05*(sample.score + expertiseCoefficient);
    }

    if( ReadyToProduce(player , sample) ){
        return 0.85*(sample.score + expertiseCoefficient);
    }

    return minScore;
};


Robot* robotForCompare = nullptr;
State* stateForCompare = nullptr;

bool compareSample(genetic_sample& s1, genetic_sample& s2){
    return sampleScore(*robotForCompare, *stateForCompare, s1) > sampleScore(*robotForCompare, *stateForCompare, s2);
}

double evalRobot(State state, Robot player, int depth) {
    if (state.turn + depth >= 200) {
        return player.score;
    }
    double evaluationScore = player.score;
    evaluationScore += expertiseCoefficient * player.totalExpertise();

    double totalSampleScore = 0;
    for (int i = 0; i < player.samples.size(); i++) {
        totalSampleScore = totalSampleScore + evalSample(state , player, player.samples[i]);
        if(!player.samples[i].diagnosed) {
            evaluationScore -= 0.01* (player.eta + Distances[player.location][DIAGNOSIS]);
        }
    }

    robotForCompare = &player;
    stateForCompare = &state;
    sort(player.samples.begin(), player.samples.end(), compareSample);  // Sorting in descending value (notice the > in compareSample)
 
    //TODO sort the samples and give 1e-2*pow(0.5,i) penalty for each missing molecule to complete sample i.
    int i = 0;
    for (vector<genetic_sample>::iterator it = player.samples.begin(); it != player.samples.end(); it++) {
        evaluationScore -= (1e-2)*pow(0.5,i)*MissingMolecules(player, *it);
        i++;
    }
    
    evaluationScore += (1e-2/4)*(10 -(player.totalMolecules()));

    return evaluationScore + totalSampleScore;
};

double eval(State state, int depth){
    return evalRobot(state , state.players[0], depth) - evalRobot(state , state.players[1], depth);
};


vector<Action> possibleMoves(State state, Robot player) {
    vector<Action> actions;
    
    if ( player.location == START ) {
        actions.push_back({GOTO, SAMPLES});
    } else if ( player.eta > 0 ) {
        actions.push_back({WAIT, 5});
    } else if( player.location == SAMPLES) {
        if ( player.samples.size() < 3) {
            if ( player.totalExpertise() < 9) {
                actions.push_back({CONNECT, 1});
            } else if ( player.totalExpertise() < 12) {
                actions.push_back({CONNECT, 2});
            } else {
                actions.push_back({CONNECT, 3});
            }
        } else {
            actions.push_back({GOTO, DIAGNOSIS});
        }
    } else if ( player.location == DIAGNOSIS ) {
        int playerSamplesCount = player.samples.size();
        // If the player has an undiagnosed sample the only move is to diagnose it, 
        //I didn't think there would be many situations where any other move would be better
        for (int i = 0; i < playerSamplesCount; i++) {
            if ( !player.samples[i].diagnosed) {
                actions.push_back({CONNECT, player.samples[i].id});
                return actions;
            }
        }
        //If the player has 3 samples, a possible move is to go to MOLECULES.
        if ( playerSamplesCount == 3 ) {
            actions.push_back({GOTO, MOLECULES});
        }

        // Pending a lot of rules here
    } else if( player.location == MOLECULES ) {
        if ( player.hasSpace() ) {
            for(int i=0;i<5;++i){
                if (state.Available[i] > 0 ) {
                    actions.push_back({CONNECT, i});
                }
            }
        }
        bool readyToProduceSamples = false;
        for(int i=0;i<player.samples.size();++i) {
            if ( ReadyToProduce(player, player.samples[i]) ) {
                readyToProduceSamples = true;
                // Does this break from for loop or if????? => break only breaks loops. What would "breaking the if" mean????
                break;
            }
        }
        if(readyToProduceSamples) {
            actions.push_back({GOTO, LABORATORY});
        } else if (player.samples.size() < 3){
            actions.push_back({GOTO, SAMPLES});
        } else {
            actions.push_back({GOTO, MOLECULES});   // Stay where you are until you have space
        }

        // Pending to do If the enemy isn't at SAMPLES and player.score + score_of_ready_samples
        // > enemy.score + score_of_enemy_produceable_samples then a possible move is to wait at the MOLECULES.

    } else if( player.location == LABORATORY) {
        bool readyToProduceSamples = false;
        int sampleToProduceId = -1;
        for(int i=0;i<player.samples.size();++i) {
            if ( ReadyToProduce(player, player.samples[i]) ) {
                readyToProduceSamples = true;
                sampleToProduceId = player.samples[i].id;
                break;
            }
        }
        if ( readyToProduceSamples ) {
            actions.push_back({CONNECT, sampleToProduceId});

            Robot enemyPlayer = state.players[!player._id];
            //TODO If a sample is ready, the enemy isn't at SAMPLES and I'm currently winning by the same criteria 
            // as in MOLECULES then a possible move is to wait at LABORATORY instead of handing the sample right away.

        } else {
            actions.push_back({GOTO, MOLECULES});
            actions.push_back({GOTO, SAMPLES});

            bool producibleSamples = false;
            for(int i=0; i < state.Samples.size();i++){
                if ( CouldProduce(player, state, state.Samples[i])){
                    producibleSamples = true;
                    break;
                }
            }
            if ( state.Samples.size() + player.samples.size() > 2){
                actions.push_back({GOTO, DIAGNOSIS});
            }
        }
        
    }

    if ( actions.size() == 0) {
        //cerr << "FUCKUP HAS HAPPENED on location: " << player.location << " OF PLAYER " << player._id <<  endl;
        actions.push_back({GOTO, SAMPLES});  // Default move when we don't know what to do
    }
    
    return actions;
};


variation miniMax(State state , int depth, int maxDepth, double alpha, double beta, Timeout& timeout) {
    if (timeout.isElapsed()) {
        //cerr << "Throwing exception";
        throw " miniMax out of time ";
    }

    if (depth == maxDepth) {
        return variation{eval(state, depth), {}};
    }
    array<vector<Action>, 2> Branch = {possibleMoves(state, state.players[0]), possibleMoves(state, state.players[1])};

    variation Best_var={-DBL_MAX,{}};
    for(Action mv:Branch[0]){
        variation Best_var2{+DBL_MAX,{}};
        double localBeta=beta;
        for(Action mv2:Branch[1]){
            State state2 = state;
            Simulate(state2, {mv, mv2});
            variation var=miniMax(state2,depth+1,maxDepth,alpha,localBeta, timeout);
            if(var.score < Best_var2.score){
                Best_var2=var;
                // This might be wrong Below
                Best_var2.Moves.insert(Best_var2.Moves.begin(),{mv,mv2});
            }
            localBeta = min(var.score , localBeta);
            if(localBeta <= alpha){
                break;
            }
        }
        if(Best_var2.score>Best_var.score){
            Best_var=Best_var2;
        }
        alpha=max(alpha,Best_var2.score);
        if(beta <= alpha){
            break;
        }
    }

    return Best_var;
};

struct Agent
{
    State state;
    Action bestAction;

    array<list<sample_template>, 3> SamplePool;

    void read();
    void think();
    void print();
};

void Agent::read()
{
    ++state.turn;
    for (int i = 0; i < 2; i++) {
        Robot& robot = state.players[i];

        string target;
        int eta, score, storageA, storageB, storageC, storageD, storageE, expertiseA, expertiseB, expertiseC, expertiseD, expertiseE;

        cin >> target >> eta >> score >> storageA >> storageB >> storageC >> storageD >> storageE >> expertiseA >> expertiseB >> expertiseC >> expertiseD >> expertiseE; cin.ignore();

        if ( target == "SAMPLES") {
            robot.location = SAMPLES;
        }
        else if ( target == "DIAGNOSIS") {
            robot.location = DIAGNOSIS;
        }
        else if ( target == "MOLECULES") {
            robot.location = MOLECULES;
        }
        else if ( target == "LABORATORY") {
            robot.location = LABORATORY;
        }
        else if ( target == "START_POS") {
            robot.location = START;
        };

        robot.eta = eta;
        robot.score = score;
        robot.Molecules[0] = storageA;
        robot.Molecules[1] = storageB;
        robot.Molecules[2] = storageC;
        robot.Molecules[3] = storageD;
        robot.Molecules[4] = storageE;
        robot.Expertise[0] = expertiseA;
        robot.Expertise[1] = expertiseB;
        robot.Expertise[2] = expertiseC;
        robot.Expertise[3] = expertiseD;
        robot.Expertise[4] = expertiseE;
        robot.samples.clear();
        robot._id = i;
    }
    int availableA , availableB, availableC, availableD, availableE;
    cin >> availableA >> availableB >> availableC >> availableD >> availableE; cin.ignore();


    state.Available[0] = availableA;
    state.Available[1] = availableB;
    state.Available[2] = availableC;
    state.Available[3] = availableD;
    state.Available[4] = availableE;

    int sampleCount;
    cin >> sampleCount; cin.ignore();
    
    state.Samples.clear();
    for (int i = 0; i < sampleCount; i++) {
        int sampleId;
        int carriedBy;
        int rank;
        string expertiseGain;
        int health;
        int costA;
        int costB;
        int costC;
        int costD;
        int costE;
        cin >> sampleId >> carriedBy >> rank >> expertiseGain >> health >> costA >> costB >> costC >> costD >> costE; cin.ignore();
        state.Samples.push_back(SamplePool[rank-1].front());
        //genetic_sample& samp = state.Samples.back();

        state.Samples[0].id = sampleId;
        state.Samples[0].owner = carriedBy;
        state.Samples[0].rank = rank;

        if(expertiseGain == "A") {
            state.Samples[0].exp = 0;
        } else if ( expertiseGain == "B") {
            state.Samples[0].exp = 1;
        } else if ( expertiseGain == "C") {
            state.Samples[0].exp = 2;
        } else if ( expertiseGain == "D") {
            state.Samples[0].exp = 3;
        } else if ( expertiseGain == "E") {
            state.Samples[0].exp = 4;
        }
        state.Samples[0].score = health;
        state.Samples[0].Cost = {costA, costB, costC, costD, costE};
        state.Samples[0].diagnosed = costA == -1 ? false : true;
        state.players[carriedBy].samples.push_back(state.Samples[0]);
    }
}


void Agent::print()
{
    
    if ( bestAction.type == GOTO ) 
    {
        cout << "GOTO " << locationToStr[bestAction.id] << endl;
    } 
    else if ( bestAction.type == ActionType::CONNECT ) 
    {
        if (state.players[0].location == MOLECULES) {
            cout << "CONNECT " << typeToStr[bestAction.id] << endl;
        } 
        else  {
            cout << "CONNECT " << bestAction.id << endl;
        }
    } 
    else if ( bestAction.type == ActionType::WAIT ) 
    {
        cout << "WAIT" << endl;
    }
}




void Agent::think()
{
    // cerr << "Current Turn: " << state.turn << endl;
    variation best_var;
    Timeout timeout(0.030);

    int depth=1;
    while(depth <= 10){
        try {
            best_var= miniMax(state, 0, depth, -DBL_MAX, +DBL_MAX, timeout);//Minimax throws an exception if time runs out
        } catch(...) {
            //cerr << "miniMax stopped at depth " << depth;
            break;
        }
        ++depth;
    }

    bestAction.type = best_var.Moves[0][0].type;
    bestAction.id = best_var.Moves[0][0].id;
    // cerr << "Max Depth Reached: " << depth << endl;
    // cerr << "BEST ACTION TYPE: " << best_var.Moves[0][0].type << endl;
    // cerr << "BEST ACTION ID: " << best_var.Moves[0][0].id << endl;
    // cerr << "Current Variation Score: " << best_var.score << endl;
}

int main()
{
    Agent agent;

    int projectCount;
    cin >> projectCount; cin.ignore();
    for (int i = 0; i < projectCount; i++) {
        int a;
        int b;
        int c;
        int d;
        int e;
        cin >> a >> b >> c >> d >> e; cin.ignore();
        agent.state.projects.push_back(Project());
        agent.state.projects[0].Target = {a, b, c, d, e};
    }


    while(1) 
    {
        agent.read();
        agent.think();
        agent.print();
    }
}
