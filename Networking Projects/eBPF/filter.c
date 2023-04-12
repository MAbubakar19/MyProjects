#define KBUILD_MODNAME "bcc"
#include <linux/bpf.h>
#include <linux/if_ether.h>
#include <linux/ip.h>
#include <linux/in.h>
#include <linux/udp.h>


int udpfilter(struct xdp_md *ctx) 
{
  int var;
  bpf_trace_printk("got a packet ");
  void *data = (void *)(long)ctx->data;
  void *data_end = (void *)(long)ctx->data_end;
  struct ethhdr *eth = data;
  if ((void*)eth + sizeof(*eth) <= data_end) 
  {
    struct iphdr *ip = data + sizeof(*eth);
    if ((void*)ip + sizeof(*ip) <= data_end) 
    {
      if (ip->protocol == IPPROTO_UDP) 
      {
        struct udphdr *udp = (void*)ip + sizeof(*ip);
        if ((void*)udp + sizeof(*udp) <= data_end) 
        {

          if (udp->dest == ntohs(20000) && (udp->source%3==0)) 
          {
            udp->dest = ntohs(20001);  
            bpf_trace_printk("YOUR PKT SEND AT 1st SERVER "); 
          }
          if ((udp->dest == ntohs(20000)) && (udp->source%3==1))
          {
            udp->dest = ntohs(20002);
            bpf_trace_printk("YOUR PKT SEND AT 2nd SERVER ");
          }
          if ((udp->dest == ntohs(20000)) && (udp->source%3==2))
          {
            udp->dest = ntohs(20003);
            bpf_trace_printk("YOUR PKT SEND AT 3rd SERVER ");
          }
        }
      }
    }
  }
  return XDP_PASS;
}
