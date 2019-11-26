/*
* This File token operations
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include "tokens.h"

Token *AddToken(BaseFloat tot_cost, BaseFloat extra_cost, ForwardLink *links, Token *next)
{
	Token *tok = new Token();

	tok->tot_cost = tot_cost;
	tok->extra_cost = extra_cost;
	tok->links = links;
	tok->next = next;

	return tok;
}

ForwardLink *AddForwardLink(Token *next_tok, int32 ilabel, int32 olabel,
				BaseFloat graph_cost, BaseFloat acoustic_cost, ForwardLink *next)
{
	struct ForwardLink * f_link = new ForwardLink();

	f_link->next_tok = next_tok;
	f_link->ilabel = ilabel;
	f_link->olabel = olabel;
	f_link->graph_cost = graph_cost;
	f_link->acoustic_cost = acoustic_cost;
	f_link->next = next;

	return f_link;
}